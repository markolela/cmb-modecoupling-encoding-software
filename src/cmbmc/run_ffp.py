"""FFP10 interpretative run, binding Section 5.

This runner:
- Validates the manifest snapshot SHA256 ledger.
- Validates per-entry file SHA256.
- Enforces manifest path safety rules.
- Applies the same binding pipeline as primary:
  load -> smooth(10') -> mask-to-0 -> patch -> encode -> aggregate -> theta/PP.
- Uses the Gaussian null reference from a primary artifact bundle for CI95
  (interpretative only, not a decision gate).

Outputs (binding names from plan):
- ffp_metrics.csv
- ffp_summary.json
- ffp_runinfo.json
- ffp_manifest_used.json
- ffp_manifest_ledger_sha256.txt
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
from typing import Any

import numpy as np
import healpy as hp

from cmbmc.smoothing import smooth_10arcmin, ELL_MAX_USED
from cmbmc.patch_centers import load_patch_centers
from cmbmc.patch_extract import extract_patch
from cmbmc.scale_quantize import coarse_grain_block_mean, quantize_u8_patch_internal
from cmbmc.gzip_proxy import gzip_smoketest_or_raise, baseline_bpc0_by_scale, bpc_from_quantized_u8, kappa_from_bpc
from cmbmc.aggregate_trend import median_aggregate_kappa, fit_theta, detrend, pp_from_kappa_theta, effects_with_ci95


S_LIST = [1, 2, 4]

MASK_PATH = Path("data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")
MASK_LEDGER = Path("preregistration/external_inputs/sha256_mask_planck_pr3_common_int_2048_r3_00.txt")

FFP_ROOT = Path("data/external/ffp10")


@dataclass(frozen=True)
class ManifestEntry:
    r: int
    local_path: str
    sha256: str
    field: int


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_single_hex(p: Path) -> str:
    s = p.read_text(encoding="utf-8").strip()
    if not s or any(c not in "0123456789abcdef" for c in s.lower()) or len(s) != 64:
        raise RuntimeError(f"Ledger file must contain exactly one SHA256 hex digest: {p.as_posix()}")
    return s.lower()


def _load_mask_binarized() -> np.ndarray:
    if not MASK_PATH.exists():
        raise FileNotFoundError(MASK_PATH.as_posix())
    if not MASK_LEDGER.exists():
        raise FileNotFoundError(MASK_LEDGER.as_posix())

    got = _sha256_file(MASK_PATH)
    exp = _read_single_hex(MASK_LEDGER)
    if got != exp:
        raise RuntimeError(f"Mask SHA256 mismatch. GOT={got} EXP={exp}")

    m = hp.read_map(MASK_PATH.as_posix(), field=0, nest=False, dtype=np.float64, verbose=False)
    m = np.asarray(m, dtype=np.float64, order="C")
    return m


def _apply_mask_to_zero(map_in: np.ndarray, mask_map: np.ndarray) -> np.ndarray:
    out = np.asarray(map_in, dtype=np.float64, order="C").copy()
    out[mask_map < np.float64(0.999)] = np.float64(0.0)
    return out


def load_ffp_manifest(manifest_path: Path) -> list[ManifestEntry]:
    obj = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries_raw = obj["entries"] if isinstance(obj, dict) and "entries" in obj else obj
    if not isinstance(entries_raw, list):
        raise RuntimeError("FFP manifest must be a list or a dict with key 'entries'.")

    out: list[ManifestEntry] = []
    for e in entries_raw:
        out.append(
            ManifestEntry(
                r=int(e["r"]),
                local_path=str(e["local_path"]),
                sha256=str(e["sha256"]).lower(),
                field=int(e.get("field", 0)),
            )
        )
    return out


def _enforce_manifest_path_safety(local_path: str) -> None:
    p = Path(local_path)
    if p.is_absolute():
        raise RuntimeError(f"Manifest local_path must not be absolute: {local_path}")
    parts = p.parts
    if ".." in parts:
        raise RuntimeError(f"Manifest local_path must not contain '..': {local_path}")
    if not local_path.startswith("data/external/ffp10/"):
        raise RuntimeError(f"Manifest local_path must be under data/external/ffp10/: {local_path}")


def _validate_manifest_snapshot(manifest_path: Path, ledger_path: Path) -> str:
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path.as_posix())
    if not ledger_path.exists():
        raise FileNotFoundError(ledger_path.as_posix())

    got = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    exp = _read_single_hex(ledger_path)
    if got != exp:
        raise RuntimeError(f"Manifest SHA256 mismatch. GOT={got} EXP={exp}")
    return got


def _load_primary_null_arrays(primary_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load theta and PP arrays from a primary artifact bundle.

    This is interpretative support for Section 5 only.
    It accepts column name variants to remain robust across artifact schema versions.
    """
    p = Path(primary_dir) / "null_effects.csv"
    if not p.exists():
        raise FileNotFoundError(p)

    theta_list: list[float] = []
    pp_list: list[float] = []

    with p.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        fields = r.fieldnames or []

        theta_key = "theta_null" if "theta_null" in fields else ("theta" if "theta" in fields else None)
        pp_key = "pp_null" if "pp_null" in fields else ("pp" if "pp" in fields else None)

        if theta_key is None or pp_key is None:
            raise KeyError(
                f"null_effects.csv has unexpected columns. "
                f"Need theta_null/theta and pp_null/pp. Got: {fields}"
            )

        for row in r:
            theta_list.append(float(row[theta_key]))
            pp_list.append(float(row[pp_key]))

    theta = np.asarray(theta_list, dtype=np.float64)
    pp = np.asarray(pp_list, dtype=np.float64)

    if theta.size == 0 or pp.size == 0:
        raise ValueError("null_effects.csv is empty or unreadable.")

    return theta, pp

def run_ffp(
    outdir: Path,
    *,
    manifest_path: Path,
    manifest_ledger_path: Path,
    primary_artifacts_dir: Path,
    baseline_seed: int = 2026022401,
) -> None:
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if any(outdir.iterdir()):
        raise RuntimeError(f"Output directory must be empty: {outdir.as_posix()}")

    t0 = time.time()

    gzip_smoketest_or_raise()

    manifest_sha256 = _validate_manifest_snapshot(manifest_path, manifest_ledger_path)
    entries = load_ffp_manifest(manifest_path)

    # Selection rule (binding): first N_FFP entries, N_FFP=100.
    entries = entries[:100]
    if len(entries) != 100:
        raise RuntimeError(f"Expected 100 entries for v2 snapshot. Got {len(entries)}")

    # Manifest snapshot rule (binding).
    (outdir / "ffp_manifest_used.json").write_bytes(manifest_path.read_bytes())
    used_sha = hashlib.sha256((outdir / "ffp_manifest_used.json").read_bytes()).hexdigest()
    (outdir / "ffp_manifest_ledger_sha256.txt").write_text(used_sha + "\n", encoding="utf-8")

    # Load null reference arrays from primary artifacts.
    theta_null, pp_null = _load_primary_null_arrays(primary_artifacts_dir)

    # Prepare fixed parts.
    bpc0 = baseline_bpc0_by_scale(baseline_seed=baseline_seed)
    centers = load_patch_centers()
    mask_map = _load_mask_binarized()

    rows: list[dict[str, Any]] = []

    for e in entries:
        _enforce_manifest_path_safety(e.local_path)
        p = Path(e.local_path)
        if not p.exists():
            raise FileNotFoundError(p.as_posix())
        got = _sha256_file(p)
        if got != e.sha256:
            raise RuntimeError(f"FFP file SHA256 mismatch for {p.as_posix()}. GOT={got} EXP={e.sha256}")

        m = hp.read_map(p.as_posix(), field=int(e.field), nest=False, dtype=np.float64, verbose=False)
        m = np.asarray(m, dtype=np.float64, order="C")

        m = smooth_10arcmin(m, ell_max_used=ELL_MAX_USED)
        m = _apply_mask_to_zero(m, mask_map)

        K = np.zeros((len(centers), len(S_LIST)), dtype=np.float64)

        for pi, c in enumerate(centers):
            A = extract_patch(m, theta0=c.theta0, phi0=c.phi0)
            for ti, s in enumerate(S_LIST):
                A_s = coarse_grain_block_mean(A, int(s))
                q = quantize_u8_patch_internal(A_s)
                bpc = bpc_from_quantized_u8(q)
                K[pi, ti] = kappa_from_bpc(bpc, bpc0[int(s)])

        ktilde = median_aggregate_kappa(K)
        _a, theta = fit_theta(ktilde)
        ktheta = detrend(ktilde, theta)
        pp = pp_from_kappa_theta(ktheta)

        rows.append(
            dict(
                r=int(e.r),
                file=p.name,
                sha256=str(e.sha256),
                theta_ffp=float(theta),
                pp_ffp=float(pp),
            )
        )

    # Write ffp_metrics.csv (per-realization).
    with (outdir / "ffp_metrics.csv").open("w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["r", "file", "sha256", "theta_ffp", "pp_ffp"])
        wr.writeheader()
        for row in rows:
            wr.writerow(row)

    # Ensemble medians (interpretative).
    theta_ffp_ref = np.median(np.asarray([r["theta_ffp"] for r in rows], dtype=np.float64))
    pp_ffp_ref = np.median(np.asarray([r["pp_ffp"] for r in rows], dtype=np.float64))

    eff = effects_with_ci95(theta_data=np.float64(theta_ffp_ref), pp_data=np.float64(pp_ffp_ref), theta_null=theta_null, pp_null=pp_null)

    summary = dict(
        n_ffp=int(len(rows)),
        s_list=S_LIST,
        baseline_seed=int(baseline_seed),
        manifest_path=str(manifest_path.as_posix()),
        manifest_sha256=str(manifest_sha256),
        manifest_used_sha256=str(used_sha),
        primary_artifacts_dir=str(primary_artifacts_dir.as_posix()),
        theta_ffp_ref=float(theta_ffp_ref),
        pp_ffp_ref=float(pp_ffp_ref),
        theta_null_ref=float(eff.theta_null_ref),
        pp_null_ref=float(eff.pp_null_ref),
        delta_theta=float(eff.delta_theta),
        ci95_delta_theta=[float(eff.ci95_delta_theta[0]), float(eff.ci95_delta_theta[1])],
        delta_pp=float(eff.delta_pp),
        ci95_delta_pp=[float(eff.ci95_delta_pp[0]), float(eff.ci95_delta_pp[1])],
    )
    (outdir / "ffp_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    runinfo = dict(
        git_commit=_git_commit(),
        platform=platform.platform(),
        python_version=_pyver(),
        numpy_version=np.__version__,
        healpy_version=hp.__version__,
        elapsed_seconds=float(time.time() - t0),
    )
    (outdir / "ffp_runinfo.json").write_text(json.dumps(runinfo, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _git_commit() -> str:
    import subprocess
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return ""


def _pyver() -> str:
    import sys
    return sys.version.replace("\n", " ")


