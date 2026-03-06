"""Primary run runner, Sections 2-4.

This module wires together the already binding-implemented components:
- Section 2.1 map loading with SHA256 validation
- Section 2.2.2 beam smoothing (10 arcmin)
- Section 2.2.1 mask binarization and mask-to-zero
- Section 2.3 patch centers (n=256, CSV order)
- Section 2.4 patch extraction (gnomonic inverse + hp.get_interp_val)
- Sections 2.5-2.6 coarse graining and quantization
- Section 2.7 gzip proxy + baseline bpc0 + smoketest
- Sections 2.8-2.10 aggregation, trend, detrending, effects and CI95

This file defines a deterministic run function and the artifact writers.
It is not executed on import.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import platform
import subprocess
import sys
import time
import zlib

import numpy as np
import healpy as hp

from cmbmc.io_planck import (
    primary_inputs_and_ledgers,
    load_map_checked,
    load_mask_checked,
    binarize_mask,
    apply_mask_to_zero,
    load_cl_tt_from_planck_theory,
    sha256_file,
)
from cmbmc.smoothing import smooth_10arcmin, ELL_MAX_USED
from cmbmc.patch_centers import load_patch_centers
from cmbmc.patch_extract import extract_patch
from cmbmc.scale_quantize import coarse_grain_block_mean, quantize_u8_patch_internal
from cmbmc.gzip_proxy import gzip_smoketest_or_raise, baseline_bpc0_by_scale, bpc_from_quantized_u8, kappa_from_bpc, S_LIST
from cmbmc.aggregate_trend import median_aggregate_kappa, fit_theta, detrend, pp_from_kappa_theta, effects_with_ci95


BASELINE_SEED = 2026022401
NULL_SEED = 2026022402
N_NULL = 512


@dataclass(frozen=True)
class PrimaryConfig:
    baseline_seed: int = BASELINE_SEED
    null_seed: int = NULL_SEED
    n_null: int = N_NULL
    ell_max_used: int = ELL_MAX_USED


def _git_commit(repo_root: Path) -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo_root.as_posix())
        return out.decode("utf-8").strip()
    except Exception:
        return ""


def _versions_snapshot() -> dict:
    import numpy
    import scipy
    import astropy
    import healpy

    return {
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "zlib_version": zlib.ZLIB_VERSION,
        "numpy_version": numpy.__version__,
        "scipy_version": scipy.__version__,
        "healpy_version": healpy.__version__,
        "astropy_version": astropy.__version__,
    }


def _ensure_new_outdir(outdir: Path) -> None:
    if outdir.exists():
        raise FileExistsError(f"Output directory already exists: {outdir.as_posix()}")
    outdir.mkdir(parents=True, exist_ok=False)


def _write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _compute_kappa_per_patch(
    *,
    map_m: np.ndarray,
    centers,
    bpc0_by_s: dict[int, np.float64],
) -> tuple[np.ndarray, list[dict]]:
    K = np.zeros((len(centers), 3), dtype=np.float64)
    per_patch_rows: list[dict] = []

    for p, c in enumerate(centers):
        A = extract_patch(map_m, theta0=c.theta0, phi0=c.phi0)

        for t, s in enumerate(S_LIST):
            A_s = coarse_grain_block_mean(A, s)
            q = quantize_u8_patch_internal(A_s)
            bpc = bpc_from_quantized_u8(q)
            bpc0 = bpc0_by_s[int(s)]
            kappa = kappa_from_bpc(bpc, bpc0)

            K[p, t] = float(kappa)

            per_patch_rows.append({
                "patch_index": int(p),
                "s": int(s),
                "bpc": float(bpc),
                "bpc0": float(bpc0),
                "kappa": float(kappa),
            })

    return K, per_patch_rows


def run_primary(outdir: str | Path, *, cfg: PrimaryConfig = PrimaryConfig()) -> None:
    """Run Sections 2-4 for both SMICA and NILC, using one shared Gaussian null ensemble.

    Artifacts written to outdir:
    metrics_per_patch.csv
    metrics.csv
    null_effects.csv
    summary.json
    null_summary.json
    runinfo.json
    """
    repo = Path(__file__).resolve().parents[2]
    out = Path(outdir).expanduser().resolve()

    _ensure_new_outdir(out)

    t0 = time.time()

    gzip_smoketest_or_raise()

    inp, led = primary_inputs_and_ledgers()

    smica = load_map_checked(map_path=inp.smica_full, sha256_ledger=led.smica_full)
    nilc = load_map_checked(map_path=inp.nilc_full, sha256_ledger=led.nilc_full)

    mask_map = load_mask_checked(mask_path=inp.mask_common, sha256_ledger=led.mask_common)
    unmasked = binarize_mask(mask_map)

    smica = smooth_10arcmin(smica, ell_max_used=cfg.ell_max_used)
    nilc = smooth_10arcmin(nilc, ell_max_used=cfg.ell_max_used)

    smica = apply_mask_to_zero(smica, unmasked)
    nilc = apply_mask_to_zero(nilc, unmasked)

    centers = load_patch_centers()

    bpc0_by_s = baseline_bpc0_by_scale(baseline_seed=int(cfg.baseline_seed))

    all_patch_rows: list[dict] = []
    metrics_rows: list[dict] = []
    summaries: dict[str, dict] = {}

    for name, m in (("SMICA", smica), ("NILC", nilc)):
        K, rows = _compute_kappa_per_patch(map_m=m, centers=centers, bpc0_by_s=bpc0_by_s)
        for r in rows:
            r["product"] = name
        all_patch_rows.extend(rows)

        kappa_tilde = median_aggregate_kappa(K)
        _a, theta_data = fit_theta(kappa_tilde)
        kappa_theta = detrend(kappa_tilde, theta_data)
        pp_data = pp_from_kappa_theta(kappa_theta)

        for t, s in enumerate(S_LIST):
            metrics_rows.append({
                "product": name,
                "s": int(s),
                "kappa_tilde": float(kappa_tilde[t]),
                "kappa_theta": float(kappa_theta[t]),
            })

        summaries[name] = {
            "theta_data": float(theta_data),
            "pp_data": float(pp_data),
        }

    cl_tt = load_cl_tt_from_planck_theory(cl_path=inp.cl_theory, sha256_ledger=led.cl_theory, ell_max_used=cfg.ell_max_used)

    null_rows: list[dict] = []
    theta_null = np.zeros(int(cfg.n_null), dtype=np.float64)
    pp_null = np.zeros(int(cfg.n_null), dtype=np.float64)

    for r in range(int(cfg.n_null)):
        seed_r = int(cfg.null_seed) + int(r)

        np.random.seed(seed_r)
        m_null = hp.synfast(
            cl_tt,
            nside=2048,
            lmax=int(cfg.ell_max_used),
            pol=False,
            new=True,
            verbose=False,
            pixwin=False,
            fwhm=0.0,
        )

        m_null = np.asarray(m_null, dtype=np.float64, order="C")
        m_null = smooth_10arcmin(m_null, ell_max_used=cfg.ell_max_used)
        m_null = apply_mask_to_zero(m_null, unmasked)

        K_null, _rows_null = _compute_kappa_per_patch(map_m=m_null, centers=centers, bpc0_by_s=bpc0_by_s)
        kappa_tilde_null = median_aggregate_kappa(K_null)
        _a_null, theta_r = fit_theta(kappa_tilde_null)
        kappa_theta_null = detrend(kappa_tilde_null, theta_r)
        pp_r = pp_from_kappa_theta(kappa_theta_null)

        theta_null[r] = float(theta_r)
        pp_null[r] = float(pp_r)

        null_rows.append({
            "r": int(r),
            "seed_r": int(seed_r),
            "theta": float(theta_r),
            "pp": float(pp_r),
        })

    for name in ("SMICA", "NILC"):
        eff = effects_with_ci95(
            theta_data=np.float64(summaries[name]["theta_data"]),
            pp_data=np.float64(summaries[name]["pp_data"]),
            theta_null=theta_null,
            pp_null=pp_null,
        )
        summaries[name].update({
            "theta_null_ref": float(eff.theta_null_ref),
            "pp_null_ref": float(eff.pp_null_ref),
            "delta_theta": float(eff.delta_theta),
            "ci95_delta_theta": [float(eff.ci95_delta_theta[0]), float(eff.ci95_delta_theta[1])],
            "delta_pp": float(eff.delta_pp),
            "ci95_delta_pp": [float(eff.ci95_delta_pp[0]), float(eff.ci95_delta_pp[1])],
        })

    _write_csv(
        out / "metrics_per_patch.csv",
        all_patch_rows,
        fieldnames=["product", "patch_index", "s", "bpc", "bpc0", "kappa"],
    )
    _write_csv(
        out / "metrics.csv",
        metrics_rows,
        fieldnames=["product", "s", "kappa_tilde", "kappa_theta"],
    )
    _write_csv(
        out / "null_effects.csv",
        null_rows,
        fieldnames=["r", "seed_r", "theta", "pp"],
    )

    _write_json(out / "summary.json", summaries)
    _write_json(out / "null_summary.json", {
        "n_null": int(cfg.n_null),
        "theta_null_ref": float(np.median(theta_null)),
        "pp_null_ref": float(np.median(pp_null)),
    })

    runinfo = {
        "git_commit": _git_commit(repo),
        "backfill_status": "none",
        "backfill_addendum": "",
        "config": {
            "baseline_seed": int(cfg.baseline_seed),
            "null_seed": int(cfg.null_seed),
            "n_null": int(cfg.n_null),
            "ell_max_used": int(cfg.ell_max_used),
            "s_list": [1, 2, 4],
        },
        "versions": _versions_snapshot(),
        "external_inputs": {
            "smica_full": {"path": inp.smica_full.as_posix(), "sha256": sha256_file(inp.smica_full)},
            "nilc_full": {"path": inp.nilc_full.as_posix(), "sha256": sha256_file(inp.nilc_full)},
            "mask_common": {"path": inp.mask_common.as_posix(), "sha256": sha256_file(inp.mask_common)},
            "cl_theory": {"path": inp.cl_theory.as_posix(), "sha256": sha256_file(inp.cl_theory)},
        },
        "timing": {
            "t_start_unix": float(t0),
            "t_end_unix": float(time.time()),
            "seconds": float(time.time() - t0),
        },
    }
    _write_json(out / "runinfo.json", runinfo)
