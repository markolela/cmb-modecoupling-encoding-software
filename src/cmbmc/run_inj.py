"""Injection suite runner, binding Section 6.

Implements:
- Injection type A: multiplicative low-ell modulation.
- Alpha grid: [0.05, 0.10, 0.15], N_inj = 128 per alpha.
- Seeding rules:
  seed_T = inj_seed + 100000*alpha_index + r
  seed_M = inj_seed + 500000 + 100000*alpha_index + r
- Low-ell spectrum for M: cl_tt_low[21:] = 0
- Normalization over unmasked set U where mask >= 0.999:
  mu = mean(M[U]), sigma = std(M[U], ddof=0), abort if sigma == 0
  M_norm = (M - mu)/sigma
  T_inj = (1 + alpha*M_norm) * T_G
- Smoothing is applied to T_inj only (after forming), then mask-to-0, then patching and encoding.

Uses Gaussian null reference from a primary artifact bundle:
paper/artifacts/primary/<RUN_ID>/null_effects.csv and null_summary.json.

Outputs (plan artifact names):
- inj_metrics.csv
- inj_summary.json
- inj_runinfo.json
- inj_null_effects.csv
- inj_null_summary.json
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import json
import platform
import subprocess
import time

import numpy as np
import healpy as hp

from cmbmc.smoothing import smooth_10arcmin, ELL_MAX_USED
from cmbmc.patch_centers import load_patch_centers
from cmbmc.patch_extract import extract_patch
from cmbmc.scale_quantize import coarse_grain_block_mean, quantize_u8_patch_internal
from cmbmc.gzip_proxy import gzip_smoketest_or_raise, baseline_bpc0_by_scale, bpc_from_quantized_u8, kappa_from_bpc, S_LIST
from cmbmc.aggregate_trend import median_aggregate_kappa, fit_theta, detrend, pp_from_kappa_theta, effects_with_ci95


# --------- Frozen external inputs (Section 2.2.1 and 4.1) ---------

MASK_PATH = Path("data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")
MASK_LEDGER = Path("preregistration/external_inputs/sha256_mask_planck_pr3_common_int_2048_r3_00.txt")

CL_PATH = Path("data/external/planck_pr3/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt")
CL_LEDGER = Path("preregistration/external_inputs/sha256_cl_theory_planck_pr3_minimum_theory_r3_01.txt")

NSIDE = 2048


@dataclass(frozen=True)
class InjConfig:
    baseline_seed: int = 2026022401
    inj_seed: int = 2026022403
    n_inj: int = 128
    alpha_list: tuple[float, float, float] = (0.05, 0.10, 0.15)


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_single_hex(path: Path) -> str:
    s = path.read_text(encoding="utf-8").strip()
    if len(s) != 64:
        raise ValueError(f"Ledger file must contain a single 64-hex digest: {path}")
    return s


def _validate_sha256(path: Path, ledger: Path) -> None:
    got = _sha256_file(path)
    exp = _read_single_hex(ledger)
    if got != exp:
        raise ValueError(f"SHA256 mismatch for {path}. GOT={got} EXP={exp}")


def _load_mask_binary() -> np.ndarray:
    """Load mask and return boolean U where mask>=0.999, binding."""
    _validate_sha256(MASK_PATH, MASK_LEDGER)
    m = hp.read_map(MASK_PATH.as_posix(), field=0, nest=False, dtype=np.float64, verbose=False)
    m = np.asarray(m, dtype=np.float64)
    return (m >= np.float64(0.999))


def _apply_mask_to_zero(map_in: np.ndarray, U: np.ndarray) -> np.ndarray:
    """Masked pixels are set to 0, binding."""
    x = np.asarray(map_in, dtype=np.float64, order="C").copy()
    x[~U] = np.float64(0.0)
    return np.asarray(x, dtype=np.float64, order="C")


def _parse_cl_tt() -> np.ndarray:
    """Parse Planck theory D_ell^TT file to C_ell^TT, binding Section 4.1."""
    _validate_sha256(CL_PATH, CL_LEDGER)

    cl_tt = np.zeros(ELL_MAX_USED + 1, dtype=np.float64)
    with CL_PATH.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            t = line.strip()
            if not t:
                continue
            c0 = t[0]
            if c0 == "#":
                continue
            if not c0.isdigit():
                continue
            parts = t.split()
            if len(parts) < 2:
                continue
            ell = int(float(parts[0]))
            dtt = float(parts[1])
            if ell < 2 or ell > ELL_MAX_USED:
                continue
            cl_tt[ell] = dtt * (2.0 * np.pi) / (ell * (ell + 1))
    cl_tt[0] = 0.0
    cl_tt[1] = 0.0
    return cl_tt


def _load_primary_null_arrays(primary_artifacts_dir: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load theta_null and pp_null arrays from primary null_effects.csv.

    Column names are accepted flexibly to avoid schema drift.
    """
    ne = primary_artifacts_dir / "null_effects.csv"
    ns = primary_artifacts_dir / "null_summary.json"
    if not ne.exists():
        raise FileNotFoundError(ne)
    if not ns.exists():
        raise FileNotFoundError(ns)

    null_summary = json.loads(ns.read_text(encoding="utf-8"))

    theta = []
    pp = []
    with ne.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        if rdr.fieldnames is None:
            raise ValueError("null_effects.csv has no header")
        fns = set([x.strip() for x in rdr.fieldnames])

        theta_key = None
        for cand in ("theta_null", "theta", "theta_r"):
            if cand in fns:
                theta_key = cand
                break
        pp_key = None
        for cand in ("pp_null", "pp", "pp_r", "PP", "pp_data"):
            if cand in fns:
                pp_key = cand
                break
        if theta_key is None or pp_key is None:
            raise ValueError(f"Unsupported null_effects.csv schema. fields={sorted(fns)}")

        for row in rdr:
            theta.append(float(row[theta_key]))
            pp.append(float(row[pp_key]))

    return np.asarray(theta, dtype=np.float64), np.asarray(pp, dtype=np.float64), null_summary


def _map_to_theta_pp(map_fullsky: np.ndarray, *, U: np.ndarray, bpc0: dict[int, np.float64]) -> tuple[np.float64, np.float64]:
    """Run binding pipeline from smoothed full-sky map to (theta, PP)."""
    # Mask-to-0 happens immediately before patch extraction, binding.
    m = _apply_mask_to_zero(map_fullsky, U)

    centers = load_patch_centers()
    K = np.zeros((len(centers), len(S_LIST)), dtype=np.float64)

    for p, c in enumerate(centers):
        A = extract_patch(m, theta0=c.theta0, phi0=c.phi0)  # float64 (256,256)
        for t, s in enumerate(S_LIST):
            A_s = coarse_grain_block_mean(A, int(s))
            q = quantize_u8_patch_internal(A_s)
            bpc = bpc_from_quantized_u8(q)
            K[p, t] = kappa_from_bpc(bpc, bpc0[int(s)])

    ktilde = median_aggregate_kappa(K)
    _a, theta = fit_theta(ktilde)
    ktheta = detrend(ktilde, theta)
    pp = pp_from_kappa_theta(ktheta)
    return np.float64(theta), np.float64(pp)


def run_inj(
    outdir: Path,
    *,
    primary_artifacts_dir: Path,
    cfg: InjConfig = InjConfig(),
) -> None:
    outdir = Path(outdir).resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise RuntimeError(f"Output directory must be new or empty: {outdir}")

    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()

    # Smoketest must pass, binding.
    gzip_smoketest_or_raise()

    # Load null reference arrays.
    theta_null, pp_null, null_summary = _load_primary_null_arrays(Path(primary_artifacts_dir).resolve())

    # Copy primary null artifacts under injection names for traceability.
    (outdir / "inj_null_effects.csv").write_bytes((Path(primary_artifacts_dir) / "null_effects.csv").read_bytes())
    (outdir / "inj_null_summary.json").write_bytes((Path(primary_artifacts_dir) / "null_summary.json").read_bytes())

    # Load mask once.
    U = _load_mask_binary()

    # Baseline bpc0 once.
    bpc0 = baseline_bpc0_by_scale(baseline_seed=int(cfg.baseline_seed))

    # Load cl_tt once.
    cl_tt = _parse_cl_tt()
    cl_tt_low = cl_tt.copy()
    cl_tt_low[21:] = 0.0

    rows = []
    alpha_list = list(cfg.alpha_list)

    for alpha_index, alpha in enumerate(alpha_list):
        for r in range(int(cfg.n_inj)):
            # 1) Draw T_G.
            seed_T = int(cfg.inj_seed) + 100000 * int(alpha_index) + int(r)
            np.random.seed(seed_T)
            T_G = hp.synfast(
                cl_tt,
                nside=NSIDE,
                lmax=ELL_MAX_USED,
                pol=False,
                new=True,
                verbose=False,
                pixwin=False,
                fwhm=0.0,
            )
            T_G = np.asarray(T_G, dtype=np.float64, order="C")

            # 2) Draw M.
            seed_M = int(cfg.inj_seed) + 500000 + 100000 * int(alpha_index) + int(r)
            np.random.seed(seed_M)
            M = hp.synfast(
                cl_tt_low,
                nside=NSIDE,
                lmax=ELL_MAX_USED,
                pol=False,
                new=True,
                verbose=False,
                pixwin=False,
                fwhm=0.0,
            )
            M = np.asarray(M, dtype=np.float64, order="C")

            # 3-4) Normalize over U.
            mu = np.mean(M[U], dtype=np.float64)
            sigma = np.std(M[U], ddof=0, dtype=np.float64)
            if float(sigma) == 0.0:
                raise RuntimeError("Injection undefined: sigma==0")
            M_norm = (M - mu) / sigma

            # 5) Form injection.
            T_inj = (np.float64(1.0) + np.float64(alpha) * M_norm) * T_G
            T_inj = np.asarray(T_inj, dtype=np.float64, order="C")

            # 6) Smooth after forming, binding.
            T_inj = smooth_10arcmin(T_inj, ell_max_used=ELL_MAX_USED)

            # 7-8) Mask-to-0, patching, encoding.
            theta, pp = _map_to_theta_pp(T_inj, U=U, bpc0=bpc0)

            rows.append(
                dict(
                    alpha=float(alpha),
                    alpha_index=int(alpha_index),
                    r=int(r),
                    seed_T=int(seed_T),
                    seed_M=int(seed_M),
                    theta_inj=float(theta),
                    pp_inj=float(pp),
                )
            )

    # Write inj_metrics.csv.
    mpath = outdir / "inj_metrics.csv"
    with mpath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # Aggregate medians per alpha.
    summary = {
        "baseline_seed": int(cfg.baseline_seed),
        "inj_seed": int(cfg.inj_seed),
        "n_inj": int(cfg.n_inj),
        "alpha_list": alpha_list,
        "primary_artifacts_dir": str(Path(primary_artifacts_dir).resolve()),
        "n_null_ref": int(len(theta_null)),
        "s_list": list(S_LIST),
        "per_alpha": {},
    }

    # Precompute null refs.
    theta_null_ref = float(np.median(theta_null))
    pp_null_ref = float(np.median(pp_null))

    for alpha in alpha_list:
        th = np.asarray([x["theta_inj"] for x in rows if x["alpha"] == float(alpha)], dtype=np.float64)
        pp = np.asarray([x["pp_inj"] for x in rows if x["alpha"] == float(alpha)], dtype=np.float64)
        th_med = float(np.median(th))
        pp_med = float(np.median(pp))

        eff = effects_with_ci95(
            theta_data=np.float64(th_med),
            pp_data=np.float64(pp_med),
            theta_null=theta_null,
            pp_null=pp_null,
        )

        summary["per_alpha"][str(alpha)] = {
            "theta_inj_med": th_med,
            "pp_inj_med": pp_med,
            "theta_null_ref": theta_null_ref,
            "pp_null_ref": pp_null_ref,
            "delta_theta_inj": float(eff.delta_theta),
            "ci95_delta_theta_inj": [float(eff.ci95_delta_theta[0]), float(eff.ci95_delta_theta[1])],
            "delta_pp_inj": float(eff.delta_pp),
            "ci95_delta_pp_inj": [float(eff.ci95_delta_pp[0]), float(eff.ci95_delta_pp[1])],
        }

    # Gate is evaluated at alpha=0.10.
    gate_alpha = "0.1"
    gate = summary["per_alpha"][gate_alpha]
    summary["gate_alpha"] = 0.10
    summary["gate_success"] = bool((gate["ci95_delta_theta_inj"][0] > 0.0) and (gate["ci95_delta_pp_inj"][0] > 0.0))

    (outdir / "inj_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    # runinfo.
    def _git(cmd: list[str]) -> str:
        try:
            return subprocess.check_output(cmd, text=True).strip()
        except Exception:
            return ""

    runinfo = {
        "git_commit": _git(["git", "rev-parse", "HEAD"]),
        "platform": platform.platform(),
        "python": sys.version if "sys" in globals() else "",
        "baseline_seed": int(cfg.baseline_seed),
        "inj_seed": int(cfg.inj_seed),
        "n_inj": int(cfg.n_inj),
        "alpha_list": alpha_list,
        "primary_artifacts_dir": str(Path(primary_artifacts_dir).resolve()),
        "mask_path": str(MASK_PATH),
        "mask_sha256": _sha256_file(MASK_PATH),
        "cl_path": str(CL_PATH),
        "cl_sha256": _sha256_file(CL_PATH),
        "elapsed_sec": float(time.time() - t0),
    }
    (outdir / "inj_runinfo.json").write_text(json.dumps(runinfo, indent=2, sort_keys=True), encoding="utf-8")


