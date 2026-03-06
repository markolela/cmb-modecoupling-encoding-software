"""HMDIFF negative control runner, binding Section 7.

Implements:
- Load HM1 and HM2 maps for SMICA and NILC.
- Construct HMDIFF = 0.5 * (HM1 - HM2).
- Apply the same binding pipeline as primary:
  smooth(10') -> mask-to-0 -> patch -> encode -> aggregate -> theta/PP.
- Compute interpretative effects vs the frozen Gaussian null reference taken from
  a primary artifact bundle (null_effects.csv).

Outputs (binding names from plan):
- hmdiff_metrics.csv
- hmdiff_summary.json
- hmdiff_runinfo.json
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time

import numpy as np
import healpy as hp

from cmbmc.smoothing import smooth_10arcmin, ELL_MAX_USED
from cmbmc.patch_centers import load_patch_centers
from cmbmc.patch_extract import extract_patch
from cmbmc.scale_quantize import coarse_grain_block_mean, quantize_u8_patch_internal
from cmbmc.gzip_proxy import (
    gzip_smoketest_or_raise,
    baseline_bpc0_by_scale,
    bpc_from_quantized_u8,
    kappa_from_bpc,
)
from cmbmc.aggregate_trend import (
    median_aggregate_kappa,
    fit_theta,
    detrend,
    pp_from_kappa_theta,
    effects_with_ci95,
)


S_LIST = [1, 2, 4]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_single_hex(p: Path) -> str:
    s = p.read_text(encoding="utf-8").strip().split()
    if not s:
        raise ValueError(f"Empty ledger file: {p.as_posix()}")
    h = s[0].strip().lower()
    if not (len(h) == 64 and all(c in "0123456789abcdef" for c in h)):
        raise ValueError(f"Ledger is not a single SHA256 hex digest: {p.as_posix()}")
    return h


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _validate_sha256_or_raise(*, file_path: Path, ledger_path: Path) -> str:
    got = _sha256_file(file_path)
    exp = _read_single_hex(ledger_path)
    if got != exp:
        raise RuntimeError(
            "SHA256 mismatch.\n"
            f"FILE:   {file_path.as_posix()}\n"
            f"LEDGER: {ledger_path.as_posix()}\n"
            f"GOT:    {got}\n"
            f"EXP:    {exp}\n"
        )
    return got


def _load_map_float64(map_path: Path) -> np.ndarray:
    m = hp.read_map(
        map_path.as_posix(),
        field=0,
        nest=False,
        dtype=np.float64,
        verbose=False,
    )
    return np.asarray(m, dtype=np.float64, order="C")


def _load_mask_binary(*, mask_path: Path) -> np.ndarray:
    mask = hp.read_map(
        mask_path.as_posix(),
        field=0,
        nest=False,
        dtype=np.float64,
        verbose=False,
    )
    mask = np.asarray(mask, dtype=np.float64, order="C")
    return (mask >= np.float64(0.999))


def _apply_mask_to_zero(m: np.ndarray, umask: np.ndarray) -> np.ndarray:
    out = np.asarray(m, dtype=np.float64, order="C").copy()
    out[~umask] = np.float64(0.0)
    return out


def _load_primary_null_arrays(primary_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Robustly load theta^(r) and PP^(r) from primary null_effects.csv.

    We accept multiple possible column spellings to prevent brittle coupling.
    """
    p = primary_dir / "null_effects.csv"
    if not p.exists():
        raise FileNotFoundError(p.as_posix())

    theta_list: list[float] = []
    pp_list: list[float] = []

    with p.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        if r.fieldnames is None:
            raise RuntimeError("null_effects.csv has no header")

        # Map lower-case -> actual name
        m = {k.strip().lower(): k for k in r.fieldnames}

        def pick(cands: list[str]) -> str:
            for c in cands:
                if c.lower() in m:
                    return m[c.lower()]
            raise KeyError(f"Missing expected columns. Have: {sorted(m.keys())}")

        theta_key = pick(["theta_null", "theta", "theta_r", "theta_value"])
        pp_key = pick(["pp_null", "pp", "pp_r", "pp_value", "pp_pct", "PP"])

        for row in r:
            theta_list.append(float(row[theta_key]))
            pp_list.append(float(row[pp_key]))

    theta = np.asarray(theta_list, dtype=np.float64)
    pp = np.asarray(pp_list, dtype=np.float64)

    if theta.size == 0 or pp.size == 0 or theta.size != pp.size:
        raise RuntimeError("Primary null arrays are empty or inconsistent")

    return theta, pp


@dataclass(frozen=True)
class HMDiffConfig:
    baseline_seed: int = 2026022401
    ell_max_used: int = ELL_MAX_USED


def _compute_theta_pp_from_map(m_fullsky: np.ndarray, *, cfg: HMDiffConfig) -> tuple[np.float64, np.float64]:
    centers = load_patch_centers()
    gzip_smoketest_or_raise()
    bpc0 = baseline_bpc0_by_scale(baseline_seed=int(cfg.baseline_seed))

    K = np.zeros((len(centers), len(S_LIST)), dtype=np.float64)

    for p_idx, c in enumerate(centers):
        A = extract_patch(m_fullsky, theta0=np.float64(c.theta0), phi0=np.float64(c.phi0))
        for t, s in enumerate(S_LIST):
            A_s = coarse_grain_block_mean(A, int(s))
            q = quantize_u8_patch_internal(A_s)
            bpc = bpc_from_quantized_u8(q)
            K[p_idx, t] = float(kappa_from_bpc(bpc, bpc0[int(s)]))

    kappa_tilde = median_aggregate_kappa(K)
    _a, theta = fit_theta(kappa_tilde)
    kappa_theta = detrend(kappa_tilde, theta)
    pp = pp_from_kappa_theta(kappa_theta)
    return np.float64(theta), np.float64(pp)


def run_hmdiff(outdir: Path, *, primary_artifacts_dir: Path | None = None, cfg: HMDiffConfig | None = None) -> None:
    if cfg is None:
        cfg = HMDiffConfig()

    outdir = Path(outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise RuntimeError(f"Output directory must be empty: {outdir.as_posix()}")
    outdir.mkdir(parents=True, exist_ok=True)

    repo = repo_root()

    # Default primary artifacts dir: latest under paper/artifacts/primary
    if primary_artifacts_dir is None:
        base = repo / "paper" / "artifacts" / "primary"
        cands = sorted([p for p in base.iterdir() if p.is_dir()])
        if not cands:
            raise RuntimeError(f"No primary artifacts found under {base.as_posix()}")
        primary_artifacts_dir = cands[-1]

    primary_artifacts_dir = Path(primary_artifacts_dir).expanduser().resolve()
    theta_null, pp_null = _load_primary_null_arrays(primary_artifacts_dir)

    # External inputs and ledgers (binding names)
    mask_path = repo / "data" / "external" / "planck_pr3" / "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"
    led_mask = repo / "preregistration" / "external_inputs" / "sha256_mask_planck_pr3_common_int_2048_r3_00.txt"

    smica_hm1 = repo / "data" / "external" / "planck_pr3" / "COM_CMB_IQU-smica_2048_R3.00_hm1.fits"
    smica_hm2 = repo / "data" / "external" / "planck_pr3" / "COM_CMB_IQU-smica_2048_R3.00_hm2.fits"
    nilc_hm1 = repo / "data" / "external" / "planck_pr3" / "COM_CMB_IQU-nilc_2048_R3.00_hm1.fits"
    nilc_hm2 = repo / "data" / "external" / "planck_pr3" / "COM_CMB_IQU-nilc_2048_R3.00_hm2.fits"

    led_smica_hm1 = repo / "preregistration" / "external_inputs" / "sha256_hm_smica_hm1.txt"
    led_smica_hm2 = repo / "preregistration" / "external_inputs" / "sha256_hm_smica_hm2.txt"
    led_nilc_hm1 = repo / "preregistration" / "external_inputs" / "sha256_hm_nilc_hm1.txt"
    led_nilc_hm2 = repo / "preregistration" / "external_inputs" / "sha256_hm_nilc_hm2.txt"

    # Validate ledgers
    _validate_sha256_or_raise(file_path=mask_path, ledger_path=led_mask)
    _validate_sha256_or_raise(file_path=smica_hm1, ledger_path=led_smica_hm1)
    _validate_sha256_or_raise(file_path=smica_hm2, ledger_path=led_smica_hm2)
    _validate_sha256_or_raise(file_path=nilc_hm1, ledger_path=led_nilc_hm1)
    _validate_sha256_or_raise(file_path=nilc_hm2, ledger_path=led_nilc_hm2)

    umask = _load_mask_binary(mask_path=mask_path)

    def do_pair(name: str, p1: Path, p2: Path) -> dict:
        HM1 = _load_map_float64(p1)
        HM2 = _load_map_float64(p2)
        hmd = np.float64(0.5) * (HM1 - HM2)
        hmd = np.asarray(hmd, dtype=np.float64, order="C")

        hmd = smooth_10arcmin(hmd, ell_max_used=int(cfg.ell_max_used))
        hmd = _apply_mask_to_zero(hmd, umask)

        theta, pp = _compute_theta_pp_from_map(hmd, cfg=cfg)
        eff = effects_with_ci95(theta_data=theta, pp_data=pp, theta_null=theta_null, pp_null=pp_null)

        return dict(
            product=name,
            theta_data=float(eff.theta_data),
            theta_null_ref=float(eff.theta_null_ref),
            delta_theta=float(eff.delta_theta),
            ci95_delta_theta=[float(eff.ci95_delta_theta[0]), float(eff.ci95_delta_theta[1])],
            pp_data=float(eff.pp_data),
            pp_null_ref=float(eff.pp_null_ref),
            delta_pp=float(eff.delta_pp),
            ci95_delta_pp=[float(eff.ci95_delta_pp[0]), float(eff.ci95_delta_pp[1])],
            gate_ci_theta_lower_gt_0=bool(float(eff.ci95_delta_theta[0]) > 0.0),
        )

    t0 = time.time()
    res_smica = do_pair("SMICA", smica_hm1, smica_hm2)
    res_nilc = do_pair("NILC", nilc_hm1, nilc_hm2)
    dt = time.time() - t0

    # Conjunctive invalidation rule (mirror Section 3 structure)
    invalidates_discovery = bool(res_smica["gate_ci_theta_lower_gt_0"] and res_nilc["gate_ci_theta_lower_gt_0"])

    # Write metrics.csv
    mpath = outdir / "hmdiff_metrics.csv"
    with mpath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "product",
                "theta_data",
                "theta_null_ref",
                "delta_theta",
                "ci95_delta_theta_low",
                "ci95_delta_theta_high",
                "pp_data",
                "pp_null_ref",
                "delta_pp",
                "ci95_delta_pp_low",
                "ci95_delta_pp_high",
                "gate_ci_theta_lower_gt_0",
            ],
        )
        w.writeheader()
        for r in (res_smica, res_nilc):
            w.writerow(
                dict(
                    product=r["product"],
                    theta_data=r["theta_data"],
                    theta_null_ref=r["theta_null_ref"],
                    delta_theta=r["delta_theta"],
                    ci95_delta_theta_low=r["ci95_delta_theta"][0],
                    ci95_delta_theta_high=r["ci95_delta_theta"][1],
                    pp_data=r["pp_data"],
                    pp_null_ref=r["pp_null_ref"],
                    delta_pp=r["delta_pp"],
                    ci95_delta_pp_low=r["ci95_delta_pp"][0],
                    ci95_delta_pp_high=r["ci95_delta_pp"][1],
                    gate_ci_theta_lower_gt_0=r["gate_ci_theta_lower_gt_0"],
                )
            )

    # Write summary.json
    spath = outdir / "hmdiff_summary.json"
    spath.write_text(
        json.dumps(
            dict(
                baseline_seed=int(cfg.baseline_seed),
                s_list=S_LIST,
                primary_artifacts_dir=primary_artifacts_dir.as_posix(),
                n_null=int(theta_null.size),
                invalidates_discovery=invalidates_discovery,
                SMICA=res_smica,
                NILC=res_nilc,
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # Write runinfo.json
    git_commit = ""
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        git_commit = ""

    rpath = outdir / "hmdiff_runinfo.json"
    rpath.write_text(
        json.dumps(
            dict(
                git_commit=git_commit,
                python_version=sys.version,
                platform=platform.platform(),
                numpy_version=np.__version__,
                healpy_version=getattr(hp, "__version__", ""),
                runtime_seconds=float(dt),
                baseline_seed=int(cfg.baseline_seed),
                ell_max_used=int(cfg.ell_max_used),
                mask_path=mask_path.as_posix(),
                smica_hm1_path=smica_hm1.as_posix(),
                smica_hm2_path=smica_hm2.as_posix(),
                nilc_hm1_path=nilc_hm1.as_posix(),
                nilc_hm2_path=nilc_hm2.as_posix(),
                primary_artifacts_dir=primary_artifacts_dir.as_posix(),
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
