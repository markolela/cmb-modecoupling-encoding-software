# scripts/make_real_patches_standalone.py
# HEALPix -> 2D-Patches (Gnomonic) aus Planck/WMAP/HI4PI, nur astropy + astropy-healpix.

from __future__ import annotations
from pathlib import Path
from typing import Tuple, List, Dict, Any, Optional
import json
import numpy as np
import healpy as hp

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS
from astropy_healpix import HEALPix

# ---------------------- SO(3) Helpers ----------------------

def _rand_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """
    Zufällige SO(3) Rotation via Quaternion. Deterministisch für gegebenes rng.
    """
    u1, u2, u3 = rng.random(3)
    q1 = np.sqrt(1.0 - u1) * np.sin(2.0 * np.pi * u2)
    q2 = np.sqrt(1.0 - u1) * np.cos(2.0 * np.pi * u2)
    q3 = np.sqrt(u1) * np.sin(2.0 * np.pi * u3)
    q4 = np.sqrt(u1) * np.cos(2.0 * np.pi * u3)
    x, y, z, w = q1, q2, q3, q4
    return np.array([
        [1 - 2*(y*y + z*z),     2*(x*y - z*w),     2*(x*z + y*w)],
        [    2*(x*y + z*w), 1 - 2*(x*x + z*z),     2*(y*z - x*w)],
        [    2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ], dtype=np.float64)

def _rotate_lonlat_rad(lon_rad: np.ndarray, lat_rad: np.ndarray, R: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotiert lon und lat in Radiant mit 3x3 Rotation R auf dem Einheitskreis.
    lon_rad und lat_rad dürfen beliebige Array Formen haben.
    """
    lon = np.asarray(lon_rad, dtype=np.float64)
    lat = np.asarray(lat_rad, dtype=np.float64)

    cl = np.cos(lat)
    v = np.stack([cl * np.cos(lon), cl * np.sin(lon), np.sin(lat)], axis=-1)

    Rt = np.asarray(R, dtype=np.float64).T
    vr = v @ Rt

    x = vr[..., 0]
    y = vr[..., 1]
    z = vr[..., 2]

    lon2 = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    lat2 = np.arcsin(np.clip(z, -1.0, 1.0))
    return lon2, lat2

# ---------------------- Performance Caches ----------------------

# Cache für das Tangentialraster. Key ist (N, fov_deg).
# Das Raster ist für alle Patches gleich und hängt nur von N und FoV ab.
_TAN_UV_CACHE: dict[tuple[int, float], tuple[np.ndarray, np.ndarray]] = {}


# ---------------------- Pfade & Defaults ----------------------

REPO = Path(__file__).resolve().parents[1]
RAW  = REPO / "data" / "raw" / "astro"
OUT  = REPO / "data" / "processed" / "astro" / "patches"
OUT.mkdir(parents=True, exist_ok=True)

# Dateinamen passend zu scripts\get_astro_data.ps1
SMICA = RAW / "planck" / "COM_CMB_IQU-smica_2048_R3.00_full.fits"     # PR3
DUST  = RAW / "planck" / "COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits"  # PR2
WMAP  = RAW / "wmap"   / "wmap_ilc_9yr_v5.fits"
HI4PI = RAW / "hi4pi"  / "HI4PI_NHI_healpix_1024.fits"
NILC = RAW / "planck" / "COM_CMB_IQU-nilc_2048_R3.00_full.fits"

# 10′-Harmonisierungs-Aliase (sm10am) – die Harmonisierung erfolgt downstream.
SMICA_HARM10 = SMICA
NILC_HARM10  = NILC

# Half-Mission (HM1/HM2), 10′-harmonisiert (sm10am)
SMICA_HM1_HARM10 = RAW / "planck" / "harmonized" / "COM_CMB_IQU-smica_2048_R3.00_hm1_sm10am.fits"
SMICA_HM2_HARM10 = RAW / "planck" / "harmonized" / "COM_CMB_IQU-smica_2048_R3.00_hm2_sm10am.fits"
NILC_HM1_HARM10  = RAW / "planck" / "harmonized" / "COM_CMB_IQU-nilc_2048_R3.00_hm1_sm10am.fits"
NILC_HM2_HARM10  = RAW / "planck" / "harmonized" / "COM_CMB_IQU-nilc_2048_R3.00_hm2_sm10am.fits"

# Patch-Parameter (Re-Run: mehr Patches, engerer FoV, strengerer Galaktik-Cut)
N_PATCHES   = 64          # vorher 16
N_PIX       = 1024        # gleich lassen
FOV_DEG     = 12.0        # vorher 20.0
LAT_CUT_DEG = 60.0        # vorher 30.0
SEEDS       = dict(smica=41, dust=42, wmap=43, hi4pi=44)

# Optionaler Run-Tag, damit ältere Ergebnisse nicht überschrieben werden
RUN_TAG = ""   # "" setzen, falls überschreiben explizit gewollt

# N3: echtes SO(3) Pre Patch Surrogat.
# Setze N3_SO3_REPS auf 0, wenn du keine N3 Datasets erzeugen willst.
N3_SO3_REPS = 0
N3_SO3_SEED0 = 900_000

WRITE_STACK = False
STACK_COMPRESS = False

def _tag(name: str) -> str:
    return f"{name}__{RUN_TAG}" if RUN_TAG else name
    


# ---------------------- FITS/HEALPix Loader ----------------------

def read_healpix_map_single(path: Path) -> Tuple[np.ndarray, int, str, str]:
    """
    Robust: liest HEALPix-Map per healpy.read_map (field=0).
    - ORDERING/NESTED wird aus dem FITS-Header gelesen.
    - Gibt (values, nside, order('ring'|'nested'), frame('galactic'|'icrs')) zurück.
    """
    if not path.exists():
        raise FileNotFoundError(path)

    # Header lesen (ext=1 bevorzugt, sonst 0)
    try:
        hdr = fits.getheader(path.as_posix(), ext=1)
    except Exception:
        hdr = fits.getheader(path.as_posix(), ext=0)

    ordering = (hdr.get("ORDERING", "RING") or "RING").strip().upper()
    nest = ordering.startswith("NEST")
    # Koordinaten-Frame
    coord = (hdr.get("COORDSYS", "GAL") or "GAL").strip().upper()
    frame = "galactic" if coord.startswith("G") else "icrs"

    # Map als 1D-Vektor laden (healpy kümmert sich um RING/NESTED)
    vals = hp.read_map(path.as_posix(), field=0, nest=nest, verbose=False)
    vals = np.asarray(vals, dtype=float).ravel()

    # Maskierte oder defekte Werte markieren.
    # hp.UNSEEN ist ein plausibler Realwert, der aber physikalisch "keine Daten" bedeutet.
    bad = (~np.isfinite(vals)) | (vals == hp.UNSEEN)
    if np.any(bad):
        vals = vals.copy()
        vals[bad] = np.nan

    npix = vals.size
    nside = int(hp.npix2nside(npix))

    order_str = "nested" if nest else "ring"
    return vals, nside, order_str, frame



def hp_from_meta(nside: int, order: str, frame: str) -> HEALPix:
    frm = Galactic() if frame == "galactic" else ICRS()
    return HEALPix(nside=nside, order=order, frame=frm)


# ---------------------- Projektion & Sampling ----------------------

def center_coord(lon_deg: float, lat_deg: float, frame: str) -> SkyCoord:
    if frame == "galactic":
        return SkyCoord(l=lon_deg * u.deg, b=lat_deg * u.deg, frame=Galactic())
    else:
        return SkyCoord(ra=lon_deg * u.deg, dec=lat_deg * u.deg, frame=ICRS())

def _get_tan_uv(N: int, fov_deg: float) -> Tuple[np.ndarray, np.ndarray]:
    key = (int(N), float(fov_deg))
    hit = _TAN_UV_CACHE.get(key)
    if hit is not None:
        return hit

    half = 0.5 * float(fov_deg)
    umax = np.tan(np.deg2rad(half))

    uu = np.linspace(-umax, +umax, int(N), dtype=np.float64)
    vv = np.linspace(-umax, +umax, int(N), dtype=np.float64)
    U, V = np.meshgrid(uu, vv)

    _TAN_UV_CACHE[key] = (U, V)
    return U, V

def _gnomonic_lonlat_grid(center: SkyCoord, N: int, fov_deg: float) -> Tuple[u.Quantity, u.Quantity]:
    """
    Inverse gnomonische Projektion.
    Wir erzeugen ein uniformes Raster im Tangentialraum, dann mappen wir es auf die Sphäre.

    Definition.
    Tangentialkoordinaten (u, v) sind Ebenenkoordinaten der TAN Projektion.
    Wir wählen u und v uniform in [-tan(half), +tan(half)].

    Orientierung.
    u ist Richtung Osten, v ist Richtung Norden, bezogen auf die Sphäre am Patch Zentrum.
    """
    ctr = center
    lon0 = float(ctr.spherical.lon.to_value(u.rad))
    lat0 = float(ctr.spherical.lat.to_value(u.rad))

    cl = np.cos(lat0)
    sl = np.sin(lat0)
    co = np.cos(lon0)
    so = np.sin(lon0)

    # Einheitsvektor zum Zentrum.
    cx = cl * co
    cy = cl * so
    cz = sl

    # Lokale Basis am Zentrum.
    # Osten, Richtung steigender lon.
    ex = -so
    ey = co
    ez = 0.0

    # Norden, Richtung steigender lat.
    nx = -co * sl
    ny = -so * sl
    nz = cl

    U, V = _get_tan_uv(N=int(N), fov_deg=float(fov_deg))

    px = cx + U * ex + V * nx
    py = cy + U * ey + V * ny
    pz = cz + U * ez + V * nz

    norm = np.sqrt(px * px + py * py + pz * pz)
    px = px / norm
    py = py / norm
    pz = pz / norm

    lon = np.arctan2(py, px)
    lat = np.arcsin(pz)

    return lon * u.rad, lat * u.rad


def gnomonic_patch(
    values: np.ndarray,
    hp_obj: HEALPix,
    center: SkyCoord,
    N: int,
    fov_deg: float,
    sampler: Any,
    so3_R_inv: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Erzeugt eine N×N Patch per echter gnomonischer TAN Projektion um 'center'
    in der nativen Frame des HEALPix Objekts.

    Optional: echtes SO(3) Pre Patch Surrogat.
    Wir rotieren die Abfragekoordinaten mit R_inv.
    Das entspricht einer global rotierten Himmelssphäre f_rot(p) = f(R_inv p).
    """
    lon, lat = _gnomonic_lonlat_grid(center, N=N, fov_deg=fov_deg)

    lon_r = lon.to_value(u.rad)
    lat_r = lat.to_value(u.rad)

    if so3_R_inv is not None:
        lon_r, lat_r = _rotate_lonlat_rad(lon_r, lat_r, so3_R_inv)

    samp = sampler(lon_r, lat_r)
    return np.asarray(samp, dtype=float).reshape(N, N)


def sample_one_center(lat_cut_deg: float, rng: np.random.Generator) -> Tuple[float, float]:
    while True:
        lon = float(rng.uniform(0.0, 360.0))
        uu = float(rng.uniform(-1.0, 1.0))
        lat = float(np.degrees(np.arcsin(uu)))
        if abs(lat) >= float(lat_cut_deg):
            return lon, lat

# ---------------------- Pipeline ----------------------

def _prepare_sampler(hp_obj: HEALPix, values: np.ndarray, interp: str) -> Any:
    """
    Erzeugt eine Sampling Funktion, die lon und lat in Radiant als ndarray erwartet.
    Entscheidet einmal pro Dataset, ob bilinear genutzt werden kann.
    """
    interp = str(interp).lower().strip()

    # Nearest Neighbor Sampler. Immer vorhanden.
    def _samp_nn(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
        lon_rad = np.mod(lon_rad, 2.0 * np.pi)
        ipix = hp_obj.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad)
        return values[ipix]

    if interp == "bilinear":
        fn = getattr(hp_obj, "interpolate_bilinear_lonlat", None)
        if fn is not None:
            try:
                # Einmaliger Signatur Test.
                test_lon = np.array([0.0], dtype=np.float64) * u.rad
                test_lat = np.array([0.0], dtype=np.float64) * u.rad
                _ = fn(test_lon, test_lat, values)

                # Bilinear Sampler mit Fallback.
                def _samp_bl(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
                    lon_rad = np.mod(lon_rad, 2.0 * np.pi)
                    try:
                        return fn(lon_rad * u.rad, lat_rad * u.rad, values)
                    except Exception:
                        return _samp_nn(lon_rad, lat_rad)

                return _samp_bl
            except Exception:
                pass

    return _samp_nn

def build_dataset(
    name: str,
    infile: Path,
    n_patches: int,
    N: int,
    fov_deg: float,
    lat_cut_deg: float,
    seed: int,
    *,
    centers_override: Optional[List[Dict[str, float]]] = None,
    so3_seed: Optional[int] = None,
    derived_from: Optional[str] = None,
) -> Dict[str, Any]:
    if not infile.exists():
        print(f"SKIP {name}: Datei fehlt -> {infile}")
        return {}

    values, nside, order, frame = read_healpix_map_single(infile)
    hp_obj = hp_from_meta(nside, order, frame)

    # Sampling Modus. bilinear ist qualitativ besser, nearest ist schneller.
    sampler = _prepare_sampler(hp_obj, values, interp="bilinear")

    so3_R_inv: Optional[np.ndarray] = None
    so3_info: Optional[Dict[str, Any]] = None
    if so3_seed is not None:
        rng_so3 = np.random.default_rng(int(so3_seed))
        R = _rand_rotation_matrix(rng_so3)
        so3_R_inv = R.T
        so3_info = {
            "seed": int(so3_seed),
            "R": R.tolist(),
        }

    ds_out = OUT / name
    ds_out.mkdir(parents=True, exist_ok=True)

    # Centers Override vorbereiten. Falls gegeben, bestimmt das die Patch-Anzahl.
    centers_override_list: Optional[List[Tuple[float, float]]] = None
    if centers_override is not None:
        centers_override_list = [
            (float(c["lon_deg"]), float(c["lat_deg"])) for c in centers_override
        ]
        n_patches = int(len(centers_override_list))

    # Optional: Stack erst nach finaler n_patches-Bestimmung anlegen.
    stack: Optional[np.ndarray] = None
    if WRITE_STACK:
        stack = np.empty((n_patches, N, N), dtype=np.float32)

    rng = np.random.default_rng(seed)

    paths: List[str] = []
    centers: List[Tuple[float, float]] = []

    print(f"== {name} :: NSIDE={nside} ORDER={order} FRAME={frame} ==")

    if centers_override_list is None:
        i = 0
        tries = 0
        max_tries = int(50 * n_patches)

        while i < n_patches:
            if tries >= max_tries:
                raise RuntimeError(
                    f"Zu viele Fehlversuche beim Patch Sampling. tries={tries}, i={i}, lat_cut={lat_cut_deg}"
                )

            L, B = sample_one_center(lat_cut_deg, rng)
            ctr = center_coord(L, B, frame)

            patch = gnomonic_patch(
                values,
                hp_obj,
                ctr,
                N=N,
                fov_deg=fov_deg,
                sampler=sampler,
                so3_R_inv=so3_R_inv,
            ).astype(np.float32)

            tries += 1

            if not np.all(np.isfinite(patch)):
                continue

            fn = ds_out / f"patch_{i:02d}.npy"
            np.save(fn, patch)
            paths.append(str(fn))
            centers.append((L, B))
            if stack is not None:
                stack[i] = patch

            print(f"  patch_{i:02d}: center=({L:.2f},{B:.2f}) deg -> {fn.name}")
            i += 1

    else:
        for i, (L, B) in enumerate(centers_override_list):
            ctr = center_coord(L, B, frame)

            patch = gnomonic_patch(
                values,
                hp_obj,
                ctr,
                N=N,
                fov_deg=fov_deg,
                sampler=sampler,
                so3_R_inv=so3_R_inv,
            ).astype(np.float32)

            if not np.all(np.isfinite(patch)):
                raise RuntimeError(
                    f"Centers override patch contains NaN. dataset={name}, i={i}, center=({L:.2f},{B:.2f})"
                )

            fn = ds_out / f"patch_{i:02d}.npy"
            np.save(fn, patch)
            paths.append(str(fn))
            centers.append((L, B))
            if stack is not None:
                stack[i] = patch

            print(f"  patch_{i:02d}: center=({L:.2f},{B:.2f}) deg -> {fn.name}")

    stack_path: Optional[Path] = None
    if WRITE_STACK:
        assert stack is not None
        stack_path = ds_out / f"{name}_stack.npz"
        if STACK_COMPRESS:
            np.savez_compressed(stack_path, patches=stack)
        else:
            np.savez(stack_path, patches=stack)

    meta = dict(
        dataset=name,
        source=str(infile),
        seed=int(seed),
        derived_from=str(derived_from) if derived_from is not None else "",
        so3_rotation=so3_info,
        nside=int(nside),
        order=order,
        frame=frame,
        projection="gnomonic_tan",
        pixel_grid="uniform_in_tangent_plane",
        N=int(N),
        fov_deg=float(fov_deg),
        lat_cut_deg=float(lat_cut_deg),
        n_patches=int(n_patches),
        centers=[{"lon_deg": float(L), "lat_deg": float(B)} for (L, B) in centers],
        patches=paths,
        stack=str(stack_path) if stack_path is not None else "",
    )
    (ds_out / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


# --- Backwards-compat alias: build(...) -> build_dataset(...)
def build(*args, **kwargs):
    """
    Kompatibilitäts-Wrapper:
    erlaubt alte Aufrufe wie build(name, infile, n_patches=..., N=..., fov=..., lat_cut=...)
    und mapped sie sicher auf build_dataset(..., fov_deg=..., lat_cut_deg=...).
    """
    # sanfte Umbenennung der Parameter
    if "fov" in kwargs and "fov_deg" not in kwargs:
        kwargs["fov_deg"] = kwargs.pop("fov")
    if "lat_cut" in kwargs and "lat_cut_deg" not in kwargs:
        kwargs["lat_cut_deg"] = kwargs.pop("lat_cut")

    # Falls jemand versehentlich beides setzt, klare Fehlermeldung
    if "fov" in kwargs or "lat_cut" in kwargs:
        raise TypeError("Bitte 'fov_deg' und 'lat_cut_deg' verwenden (nicht 'fov'/'lat_cut').")

    return build_dataset(*args, **kwargs)

def build_so3_surrogate_from_base(base_dataset: str, out_dataset: str, so3_seed: int) -> Dict[str, Any]:
    """
    Baut ein echtes SO(3) Pre Patch Surrogat Dataset.
    Grundlage sind die Centers aus base_dataset.
    Die globale Rotation wird im Sampling angewendet, bevor Patch Werte gezogen werden.
    """
    base_meta_path = OUT / base_dataset / "meta.json"
    if not base_meta_path.exists():
        raise FileNotFoundError(f"Base meta fehlt: {base_meta_path}")

    base_meta = json.loads(base_meta_path.read_text(encoding="utf-8"))

    centers = base_meta.get("centers", [])
    if not centers:
        raise ValueError(f"Base meta hat keine centers: {base_meta_path}")

    src = base_meta.get("source", "")
    if not src:
        raise ValueError(f"Base meta hat kein source: {base_meta_path}")

    return build_dataset(
        name=out_dataset,
        infile=Path(src),
        n_patches=int(base_meta.get("n_patches", len(centers))),
        N=int(base_meta["N"]),
        fov_deg=float(base_meta["fov_deg"]),
        lat_cut_deg=float(base_meta["lat_cut_deg"]),
        seed=int(base_meta.get("seed", 0)),
        centers_override=centers,
        so3_seed=int(so3_seed),
        derived_from=str(base_dataset),
    )

def main():
    print("=== make_real_patches_standalone: START ===")
    print(f"Repo: {REPO}")
    print(f"RAW:  {RAW}")
    print(f"OUT:  {OUT}")

    results: Dict[str, Any] = {}

    """results["planck_smica_cmb"] = build_dataset(
        _tag("planck_smica_cmb"), SMICA,
        n_patches=N_PATCHES, N=N_PIX, fov_deg=FOV_DEG, lat_cut_deg=LAT_CUT_DEG, seed=SEEDS["smica"]
    ) or {}"""

    # --- strengere Varianten (FOV 10°, |b|>70°, 64 Patches) ---
    """results["planck_smica_cmb__fov10_b70_n64"] = build_dataset(
        "planck_smica_cmb__fov10_b70_n64", SMICA,
        n_patches=64, N=1024, fov_deg=10.0, lat_cut_deg=70.0, seed=SEEDS["smica"]
    ) or {}

    results["wmap_ilc9__fov10_b70_n64"] = build_dataset(
        "wmap_ilc9__fov10_b70_n64", WMAP,
        n_patches=64, N=1024, fov_deg=10.0, lat_cut_deg=70.0, seed=SEEDS["wmap"]
    ) or {}

    results["planck_gnilc_dust__fov10_b70_n64"] = build_dataset(
        "planck_gnilc_dust__fov10_b70_n64", DUST,
        n_patches=64, N=1024, fov_deg=10.0, lat_cut_deg=70.0, seed=SEEDS["dust"]
    ) or {}

    results["planck_gnilc_dust"] = build_dataset(
        _tag("planck_gnilc_dust"), DUST,
        n_patches=N_PATCHES, N=N_PIX, fov_deg=FOV_DEG, lat_cut_deg=LAT_CUT_DEG, seed=SEEDS["dust"]
    ) or {}

    results["planck_smica_cmb__fov12_b60_n128"] = build("planck_smica_cmb__fov12_b60_n128", SMICA, n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=41)

    results["planck_smica_cmb__fov12_b60_n64_seed42"] = build(
        "planck_smica_cmb__fov12_b60_n64_seed42",
        SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=60.0, seed=42
    )

    results["planck_smica_cmb__fov12_b60_n64_seed43"] = build(
        "planck_smica_cmb__fov12_b60_n64_seed43",
        SMICA, n_patches=64, N=1024, fov=12.0, lat_cut=60.0, seed=43
    )

    # --- SMICA fov12 |b|>60, zusätzliche Seeds ---
    results["planck_smica_cmb__fov12_b60_n64_seed44"] = build(
        "planck_smica_cmb__fov12_b60_n64_seed44",
        SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=60.0, seed=44
    )

    results["planck_smica_cmb__fov12_b60_n64_seed45"] = build(
        "planck_smica_cmb__fov12_b60_n64_seed45",
        SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=60.0, seed=45
    )

    # --- SMICA fov10 |b|>70, zusätzliche Seeds ---
    results["planck_smica_cmb__fov10_b70_n64_seed44"] = build_dataset(
        "planck_smica_cmb__fov10_b70_n64_seed44", SMICA,
        n_patches=64, N=1024, fov_deg=10.0, lat_cut_deg=70.0, seed=44
    )

    results["planck_smica_cmb__fov10_b70_n64_seed45"] = build_dataset(
        "planck_smica_cmb__fov10_b70_n64_seed45", SMICA,
        n_patches=64, N=1024, fov_deg=10.0, lat_cut_deg=70.0, seed=45
    )

    # --- SMICA: gleicher FoV (12°), strengere Galaktik-Cuts ---
    results["planck_smica_cmb__fov12_b70_n64"] = build(
        "planck_smica_cmb__fov12_b70_n64", SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=70.0, seed=SEEDS["smica"]
    )

    results["planck_smica_cmb__fov12_b75_n64"] = build(
        "planck_smica_cmb__fov12_b75_n64", SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=75.0, seed=SEEDS["smica"]
    )

    # WMAP ILC9 |b|>60°, FoV 12°, n=128 (Seeds 41–45)
    results["wmap_ilc9__fov12_b60_n128"] = build_dataset(
        "wmap_ilc9__fov12_b60_n128", WMAP,
        n_patches=128, N=1024, fov_deg=12.0, lat_cut_deg=60.0, seed=41
    )
    for sd in [42,43,44,45]:
        results[f"wmap_ilc9__fov12_b60_n128_seed{sd}"] = build_dataset(
            f"wmap_ilc9__fov12_b60_n128_seed{sd}", WMAP,
            n_patches=128, N=1024, fov_deg=12.0, lat_cut_deg=60.0, seed=sd
        )

    # SMICA |b|>65°, FoV 12°, n=128
    results["planck_smica_cmb__fov12_b65_n128"] = build(
        "planck_smica_cmb__fov12_b65_n128", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    for sd in [42,43,44,45]:
        results[f"planck_smica_cmb__fov12_b65_n128_seed{sd}"] = build(
            f"planck_smica_cmb__fov12_b65_n128_seed{sd}", SMICA,
            n_patches=128, N=1024, fov=12.0, lat_cut=65.0, seed=sd
        )

    # scripts/make_real_patches_standalone.py – Ergänzen
    results["planck_smica_cmb__fov12_b60_n256"] = build(
        "planck_smica_cmb__fov12_b60_n256", SMICA,
        n_patches=256, N=1024, fov=12.0, lat_cut=60.0, seed=41
    )
    for sd in [42,43,44,45]:
        results[f"planck_smica_cmb__fov12_b60_n256_seed{sd}"] = build(
            f"planck_smica_cmb__fov12_b60_n256_seed{sd}", SMICA,
            n_patches=256, N=1024, fov=12.0, lat_cut=60.0, seed=sd
        )

    # scripts/make_real_patches_standalone.py (zusätzliche Einträge)
    results["planck_smica_cmb__fov12_b70_n128"] = build(
        "planck_smica_cmb__fov12_b70_n128", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=70.0, seed=SEEDS["smica"]
    )
    results["planck_smica_cmb__fov12_b75_n128"] = build(
        "planck_smica_cmb__fov12_b75_n128", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=75.0, seed=SEEDS["smica"]
    )

    results["planck_smica_cmb__fov12_b60_n64"] = build(
        "planck_smica_cmb__fov12_b60_n64", SMICA,
        n_patches=64, N=1024, fov=12.0, lat_cut=60.0, seed=41
    )

    # --- SMICA |b|>60°, FoV 12°, n=128: zusätzliche Seeds ---
    results["planck_smica_cmb__fov12_b60_n128_seed42"] = build(
        "planck_smica_cmb__fov12_b60_n128_seed42", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=42
    )
    results["planck_smica_cmb__fov12_b60_n128_seed43"] = build(
        "planck_smica_cmb__fov12_b60_n128_seed43", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=43
    )
    results["planck_smica_cmb__fov12_b60_n128_seed44"] = build(
        "planck_smica_cmb__fov12_b60_n128_seed44", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=44
    )
    results["planck_smica_cmb__fov12_b60_n128_seed45"] = build(
        "planck_smica_cmb__fov12_b60_n128_seed45", SMICA,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=45
    )"""
    """
    # --- NILC sm10am: FoV 12°, |b|≥{65,70,75}, n=256, Seed=41 ---
    # Einheitliche Tags: planck_nilc_cmb__fov12_b{65,70,75}_n256_sm10am
    results["planck_nilc_cmb__fov12_b65_n256_sm10am"] = build(
        "planck_nilc_cmb__fov12_b65_n256_sm10am", NILC_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_nilc_cmb__fov12_b70_n256_sm10am"] = build(
        "planck_nilc_cmb__fov12_b70_n256_sm10am", NILC_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_nilc_cmb__fov12_b75_n256_sm10am"] = build(
        "planck_nilc_cmb__fov12_b75_n256_sm10am", NILC_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )"""
    """
    # --- NILC: FoV 12°, |b|>60°, n=128 (Seeds 41–45, analog zu SMICA) ---
    results["planck_nilc_cmb__fov12_b60_n128"] = build(
        "planck_nilc_cmb__fov12_b60_n128", NILC,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=41
    )
    results["planck_nilc_cmb__fov12_b60_n128_seed42"] = build(
        "planck_nilc_cmb__fov12_b60_n128_seed42", NILC,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=42
    )
    results["planck_nilc_cmb__fov12_b60_n128_seed43"] = build(
        "planck_nilc_cmb__fov12_b60_n128_seed43", NILC,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=43
    )
    results["planck_nilc_cmb__fov12_b60_n128_seed44"] = build(
        "planck_nilc_cmb__fov12_b60_n128_seed44", NILC,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=44
    )
    results["planck_nilc_cmb__fov12_b60_n128_seed45"] = build(
        "planck_nilc_cmb__fov12_b60_n128_seed45", NILC,
        n_patches=128, N=1024, fov=12.0, lat_cut=60.0, seed=45
    )"""
    """
    results["wmap_ilc9__fov12_b60_n64"] = build_dataset(
        "wmap_ilc9__fov12_b60_n64", WMAP,
        n_patches=64, N=1024, fov_deg=12.0, lat_cut_deg=60.0, seed=SEEDS["wmap"]
    )

    results["hi4pi_nhi"] = build_dataset(
        _tag("hi4pi_nhi"), HI4PI,
        n_patches=N_PATCHES, N=N_PIX, fov_deg=FOV_DEG, lat_cut_deg=LAT_CUT_DEG, seed=SEEDS["hi4pi"]
    ) or {}"""
    """# --- SMICA sm10am: FoV 12°, |b|≥{65,70,75}, n=256, Seed=41 ---
    # Einheitliche Tags: planck_smica_cmb__fov12_b{65,70,75}_n256_sm10am
    results["planck_smica_cmb__fov12_b65_n256_sm10am"] = build(
        "planck_smica_cmb__fov12_b65_n256_sm10am", SMICA_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_smica_cmb__fov12_b70_n256_sm10am"] = build(
        "planck_smica_cmb__fov12_b70_n256_sm10am", SMICA_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_smica_cmb__fov12_b75_n256_sm10am"] = build(
        "planck_smica_cmb__fov12_b75_n256_sm10am", SMICA_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )"""

    # Einheitliche Tags: <map>_hm{1,2}__fov12_b{65,70,75}_n256_sm10am  (Seed=41)

    # -- SMICA HM1/HM2 -----------------------------------------------------------
    results["planck_smica_cmb_hm1__fov12_b65_n256_sm10am"] = build(
        "planck_smica_cmb_hm1__fov12_b65_n256_sm10am", SMICA_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_smica_cmb_hm1__fov12_b70_n256_sm10am"] = build(
        "planck_smica_cmb_hm1__fov12_b70_n256_sm10am", SMICA_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_smica_cmb_hm1__fov12_b75_n256_sm10am"] = build(
        "planck_smica_cmb_hm1__fov12_b75_n256_sm10am", SMICA_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )

    results["planck_smica_cmb_hm2__fov12_b65_n256_sm10am"] = build(
        "planck_smica_cmb_hm2__fov12_b65_n256_sm10am", SMICA_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_smica_cmb_hm2__fov12_b70_n256_sm10am"] = build(
        "planck_smica_cmb_hm2__fov12_b70_n256_sm10am", SMICA_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_smica_cmb_hm2__fov12_b75_n256_sm10am"] = build(
        "planck_smica_cmb_hm2__fov12_b75_n256_sm10am", SMICA_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )

    # -- NILC HM1/HM2 ------------------------------------------------------------
    results["planck_nilc_cmb_hm1__fov12_b65_n256_sm10am"] = build(
        "planck_nilc_cmb_hm1__fov12_b65_n256_sm10am", NILC_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_nilc_cmb_hm1__fov12_b70_n256_sm10am"] = build(
        "planck_nilc_cmb_hm1__fov12_b70_n256_sm10am", NILC_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_nilc_cmb_hm1__fov12_b75_n256_sm10am"] = build(
        "planck_nilc_cmb_hm1__fov12_b75_n256_sm10am", NILC_HM1_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )

    results["planck_nilc_cmb_hm2__fov12_b65_n256_sm10am"] = build(
        "planck_nilc_cmb_hm2__fov12_b65_n256_sm10am", NILC_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=65.0, seed=41
    )
    results["planck_nilc_cmb_hm2__fov12_b70_n256_sm10am"] = build(
        "planck_nilc_cmb_hm2__fov12_b70_n256_sm10am", NILC_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=70.0, seed=41
    )
    results["planck_nilc_cmb_hm2__fov12_b75_n256_sm10am"] = build(
        "planck_nilc_cmb_hm2__fov12_b75_n256_sm10am", NILC_HM2_HARM10,
        n_patches=256, N=1024, fov=12.0, lat_cut=75.0, seed=41
    )

    base_tags: List[str] = []

    if int(N3_SO3_REPS) > 0:
        base_tags = [k for k, v in results.items() if isinstance(v, dict) and v]

    for base_tag in base_tags:
        base_meta_path = OUT / base_tag / "meta.json"
        if not base_meta_path.exists():
            print(f"SKIP N3: Base meta fehlt -> {base_meta_path}")
            continue

        for r in range(int(N3_SO3_REPS)):
            out_tag = f"{base_tag}__n3so3_r{r:03d}"
            results[out_tag] = build_so3_surrogate_from_base(
                base_dataset=base_tag,
                out_dataset=out_tag,
                so3_seed=int(N3_SO3_SEED0) + int(r),
            ) or {}

    manifest = OUT / "patches_manifest.json"
    manifest.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Manifest: {manifest}")
    print("=== DONE ===")


if __name__ == "__main__":
    main()
