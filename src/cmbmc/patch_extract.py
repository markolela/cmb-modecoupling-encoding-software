"""Patch extraction, binding Section 2.4.

Implements the frozen gnomonic inverse projection and sampling rule.

Binding constants.
- Projection: gnomonic
- FOV: 12 deg
- Patch grid: 256 x 256
- Array convention: A[j,i], j north (y), i east (x)
- Sampling: hp.get_interp_val with theta/phi raveled in C order, nest=False
- Output: float64, C order
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import healpy as hp


N = 256
FOV_RAD = np.deg2rad(12.0)
TWO_PI = np.float64(2.0 * np.pi)


@dataclass(frozen=True)
class PatchGrid:
    X: np.ndarray
    Y: np.ndarray
    rho: np.ndarray
    c: np.ndarray
    z: np.ndarray
    rho_safe: np.ndarray


_GRID_CACHE: PatchGrid | None = None


def _get_grid() -> PatchGrid:
    global _GRID_CACHE
    if _GRID_CACHE is not None:
        return _GRID_CACHE

    n = int(N)
    fov_rad = np.float64(FOV_RAD)
    delta = fov_rad / np.float64(n)

    i = np.arange(n, dtype=np.float64)
    j = np.arange(n, dtype=np.float64)

    half = np.float64(n) / np.float64(2.0)
    x = (i + np.float64(0.5) - half) * delta
    y = (j + np.float64(0.5) - half) * delta

    X, Y = np.meshgrid(x, y, indexing="xy")  # X[j,i]=x[i], Y[j,i]=y[j]

    rho = np.sqrt(X * X + Y * Y)
    c = np.arctan(rho)

    z = np.isclose(rho, np.float64(0.0), atol=np.float64(1e-12))
    rho_safe = np.where(z, np.float64(1.0), rho)

    _GRID_CACHE = PatchGrid(
        X=np.asarray(X, dtype=np.float64, order="C"),
        Y=np.asarray(Y, dtype=np.float64, order="C"),
        rho=np.asarray(rho, dtype=np.float64, order="C"),
        c=np.asarray(c, dtype=np.float64, order="C"),
        z=np.asarray(z, dtype=bool),
        rho_safe=np.asarray(rho_safe, dtype=np.float64, order="C"),
    )
    return _GRID_CACHE


def extract_patch(m: np.ndarray, *, theta0: np.float64, phi0: np.float64) -> np.ndarray:
    """Extract one (256,256) gnomonic patch from a full-sky HEALPix map.

    m must be a float64 array in RING ordering.
    theta0, phi0 are HEALPix angles in radians (theta colatitude, phi longitude).
    """
    grid = _get_grid()

    theta0 = np.float64(theta0)
    phi0 = np.float64(phi0)
    phi0 = np.mod(phi0, TWO_PI)

    X = grid.X
    Y = grid.Y
    c = grid.c
    z = grid.z
    rho_safe = grid.rho_safe

    cos_theta = np.cos(c) * np.cos(theta0) + (Y * np.sin(c) * np.sin(theta0)) / rho_safe
    cos_theta = np.where(z, np.cos(theta0), cos_theta)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    theta = np.arccos(cos_theta)

    denom = rho_safe * np.sin(theta0) * np.cos(c) - Y * np.cos(theta0) * np.sin(c)
    phi = phi0 + np.arctan2(X * np.sin(c), denom)

    theta = np.where(z, theta0, theta)
    phi = np.where(z, phi0, phi)

    phi = np.mod(phi, TWO_PI)

    mm = np.asarray(m, dtype=np.float64)

    vals = hp.get_interp_val(
        mm,
        theta.ravel(order="C"),
        phi.ravel(order="C"),
        nest=False,
    )

    A = np.asarray(vals.reshape((int(N), int(N)), order="C"), dtype=np.float64, order="C")
    return A
