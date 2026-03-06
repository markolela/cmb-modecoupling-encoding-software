"""Patch centers loading, binding Section 2.3.

- Loads CSV via numpy.genfromtxt only (no pandas).
- Requires exactly N_patches = 256 rows.
- Requires columns theta_rad and phi_rad.
- Enforces phi normalization: phi = mod(phi, 2*pi).
- Patch index order is CSV row order with no reordering.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


N_PATCHES = 256


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class PatchCenter:
    theta0: np.float64
    phi0: np.float64


def load_patch_centers(csv_path: str | Path = "preregistration/patch_centers/patch_centers_n256_fov12_galactic_v1.csv") -> list[PatchCenter]:
    p = Path(csv_path)
    if not p.is_absolute():
        p = repo_root() / p

    centers = np.genfromtxt(
        p.as_posix(),
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )

    # Binding validity checks
    if getattr(centers, "ndim", None) != 1:
        raise ValueError(f"Invalid centers: ndim != 1. Got ndim={getattr(centers, 'ndim', None)}")
    if len(centers) != N_PATCHES:
        raise ValueError(f"Invalid centers: expected {N_PATCHES} rows, got {len(centers)}")

    # Required columns
    names = set(centers.dtype.names or ())
    if "theta_rad" not in names or "phi_rad" not in names:
        raise ValueError(f"Invalid centers: required columns theta_rad, phi_rad. Found {sorted(names)}")

    out: list[PatchCenter] = []
    two_pi = np.float64(2.0 * np.pi)

    for k in range(N_PATCHES):
        theta0 = np.float64(centers["theta_rad"][k])
        phi0 = np.float64(centers["phi_rad"][k])
        phi0 = np.mod(phi0, two_pi)
        out.append(PatchCenter(theta0=theta0, phi0=phi0))

    return out
