"""Planck PR3 I/O and SHA256 validation.

This module implements only frozen protocol elements from the preregistration plan v3.3.
It validates external inputs via committed SHA256 ledgers and loads maps with the binding
healpy.read_map call signature.

All arrays are float64. All HEALPix maps are loaded with nest=False, field=0, verbose=False.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import re
from typing import Iterable

import numpy as np
import healpy as hp


ELL_MAX_USED = 6143


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_singleline_hex(path: Path) -> str:
    s = path.read_text(encoding="utf-8").strip()
    s = s.replace("\r", "").replace(" ", "")
    if not re.fullmatch(r"[0-9a-f]{64}", s):
        raise ValueError(f"Ledger is not a single lowercase sha256 hex digest, path={path.as_posix()}, got={s!r}")
    return s


def sha256_file(path: Path, chunk_bytes: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_bytes), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_sha256(file_path: Path, ledger_path: Path) -> None:
    if not file_path.exists():
        raise FileNotFoundError(file_path.as_posix())
    if not ledger_path.exists():
        raise FileNotFoundError(ledger_path.as_posix())
    got = sha256_file(file_path)
    exp = _read_singleline_hex(ledger_path)
    if got != exp:
        raise ValueError(
            "SHA256 mismatch. "
            f"file={file_path.as_posix()} got={got} exp={exp} ledger={ledger_path.as_posix()}"
        )


def load_map_checked(*, map_path: Path, sha256_ledger: Path) -> np.ndarray:
    """Binding disk-load call for temperature maps."""
    assert_sha256(map_path, sha256_ledger)
    m = hp.read_map(
        map_path.as_posix(),
        field=0,
        nest=False,
        dtype=np.float64,
        verbose=False,
    )
    return np.asarray(m, dtype=np.float64, order="C")


def load_mask_checked(*, mask_path: Path, sha256_ledger: Path) -> np.ndarray:
    """Binding disk-load call for the common intensity mask."""
    assert_sha256(mask_path, sha256_ledger)
    m = hp.read_map(
        mask_path.as_posix(),
        field=0,
        nest=False,
        dtype=np.float64,
        verbose=False,
    )
    return np.asarray(m, dtype=np.float64, order="C")


def binarize_mask(mask_map: np.ndarray) -> np.ndarray:
    """Binding binarization: unmasked iff mask_value >= 0.999."""
    m = np.asarray(mask_map, dtype=np.float64)
    return (m >= np.float64(0.999))


def apply_mask_to_zero(map_in: np.ndarray, unmasked: np.ndarray) -> np.ndarray:
    """Binding mask application: masked pixels set to 0. Returns float64 C-order."""
    x = np.asarray(map_in, dtype=np.float64).copy(order="C")
    u = np.asarray(unmasked, dtype=bool)
    if x.shape != u.shape:
        raise ValueError(f"Shape mismatch. map={x.shape} unmasked={u.shape}")
    x[~u] = np.float64(0.0)
    return np.asarray(x, dtype=np.float64, order="C")


def load_cl_tt_from_planck_theory(*, cl_path: Path, sha256_ledger: Path, ell_max_used: int = ELL_MAX_USED) -> np.ndarray:
    """Section 4.1. Parse Planck theory D_ell^TT to C_ell^TT.

    Binding rules.
    Read ASCII, line by line.
    Ignore empty lines.
    Ignore lines where first non-space char is not a digit.
    Ignore lines where first non-space char is '#'.
    Parse tokens as floats.
    tokens[0] is ell, tokens[1] is D_ell^TT in microK^2.
    Convert: C_ell^TT = D_ell^TT * 2*pi / (ell*(ell+1)).
    Fill cl_tt[ell], skip ell < 2, skip ell > ell_max_used.
    Set cl_tt[0]=0 and cl_tt[1]=0.
    """
    assert_sha256(cl_path, sha256_ledger)

    cl_tt = np.zeros(int(ell_max_used) + 1, dtype=np.float64)

    with cl_path.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            c0 = s[0]
            if c0 == "#":
                continue
            if not c0.isdigit():
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                ell = int(float(parts[0]))
                d_tt = float(parts[1])
            except Exception:
                continue
            if ell < 2 or ell > int(ell_max_used):
                continue
            cl_tt[ell] = d_tt * (2.0 * np.pi) / (ell * (ell + 1))

    cl_tt[0] = np.float64(0.0)
    cl_tt[1] = np.float64(0.0)
    return cl_tt


@dataclass(frozen=True)
class PrimaryInputs:
    smica_full: Path
    nilc_full: Path
    smica_hm1: Path
    smica_hm2: Path
    nilc_hm1: Path
    nilc_hm2: Path
    mask_common: Path
    cl_theory: Path


@dataclass(frozen=True)
class PrimaryLedgers:
    smica_full: Path
    nilc_full: Path
    smica_hm1: Path
    smica_hm2: Path
    nilc_hm1: Path
    nilc_hm2: Path
    mask_common: Path
    cl_theory: Path


def primary_inputs_and_ledgers() -> tuple[PrimaryInputs, PrimaryLedgers]:
    r = repo_root()

    inp = PrimaryInputs(
        smica_full=r / "data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_full.fits",
        nilc_full=r / "data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_full.fits",
        smica_hm1=r / "data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_hm1.fits",
        smica_hm2=r / "data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_hm2.fits",
        nilc_hm1=r / "data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_hm1.fits",
        nilc_hm2=r / "data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_hm2.fits",
        mask_common=r / "data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits",
        cl_theory=r / "data/external/planck_pr3/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",
    )

    led = PrimaryLedgers(
        smica_full=r / "preregistration/external_inputs/sha256_data_smica_pr3_full.txt",
        nilc_full=r / "preregistration/external_inputs/sha256_data_nilc_pr3_full.txt",
        smica_hm1=r / "preregistration/external_inputs/sha256_hm_smica_hm1.txt",
        smica_hm2=r / "preregistration/external_inputs/sha256_hm_smica_hm2.txt",
        nilc_hm1=r / "preregistration/external_inputs/sha256_hm_nilc_hm1.txt",
        nilc_hm2=r / "preregistration/external_inputs/sha256_hm_nilc_hm2.txt",
        mask_common=r / "preregistration/external_inputs/sha256_mask_planck_pr3_common_int_2048_r3_00.txt",
        cl_theory=r / "preregistration/external_inputs/sha256_cl_theory_planck_pr3_minimum_theory_r3_01.txt",
    )
    return inp, led


def validate_primary_inputs() -> None:
    inp, led = primary_inputs_and_ledgers()

    pairs: Iterable[tuple[Path, Path]] = (
        (inp.smica_full, led.smica_full),
        (inp.nilc_full, led.nilc_full),
        (inp.smica_hm1, led.smica_hm1),
        (inp.smica_hm2, led.smica_hm2),
        (inp.nilc_hm1, led.nilc_hm1),
        (inp.nilc_hm2, led.nilc_hm2),
        (inp.mask_common, led.mask_common),
        (inp.cl_theory, led.cl_theory),
    )

    for fp, lp in pairs:
        assert_sha256(fp, lp)


