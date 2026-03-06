"""Beam smoothing, binding Section 2.2.2.

This module implements the frozen healpy.smoothing call with:
FWHM = 10 arcmin, iter=3, lmax=ell_max_used, verbose=False, nest=False.
All arrays are float64 and returned as C-order contiguous arrays.
"""

from __future__ import annotations

import numpy as np
import healpy as hp

ELL_MAX_USED = 6143

# Binding definition.
FWHM_10ARCMIN_RAD = np.deg2rad(10.0 / 60.0)


def smooth_10arcmin(map_in: np.ndarray, *, ell_max_used: int = ELL_MAX_USED) -> np.ndarray:
    """Apply the binding smoothing call.

    hp.smoothing(
      map_in,
      fwhm=FWHM_10ARCMIN_RAD,
      sigma=None,
      beam_window=None,
      pol=False,
      iter=3,
      lmax=ell_max_used,
      mmax=None,
      use_weights=False,
      use_pixel_weights=False,
      datapath=None,
      verbose=False,
      nest=False
    )
    """
    m = np.asarray(map_in, dtype=np.float64, order="C")

    out = hp.smoothing(
        m,
        fwhm=float(FWHM_10ARCMIN_RAD),
        sigma=None,
        beam_window=None,
        pol=False,
        iter=3,
        lmax=int(ell_max_used),
        mmax=None,
        use_weights=False,
        use_pixel_weights=False,
        datapath=None,
        verbose=False,
        nest=False,
    )
    return np.asarray(out, dtype=np.float64, order="C")
