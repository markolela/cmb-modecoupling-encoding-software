"""Scale workflow and quantization, binding Sections 2.5 and 2.6.

Section 2.5 (coarse graining).
- Fixed scale set for this project: s_list = [1, 2, 4]
- If 256 mod s != 0: apply periodic wrap padding to make dimensions divisible by s.
- Apply s x s block-mean with fixed segmentation starting at A[0,0].
- Axis convention: axis 0 is y (rows), axis 1 is x (columns).
- Block mean:
  H = 256 // s
  W = 256 // s
  A_s = A.reshape(H, s, W, s).mean(axis=(1,3))

Section 2.6 (quantization).
- Patch-internal min/max scaling.
- If xmax == xmin: all zeros.
- Else:
  q_raw = floor(255 * (x - xmin) / (xmax - xmin))
  q_raw = clip(q_raw, 0, 255)
  q = q_raw.astype(uint8, copy=False)
- Byte order: q.ravel(order="C").tobytes()
"""

from __future__ import annotations

import numpy as np


N_PATCH = 256


def _wrap_pad_to_multiple(a: np.ndarray, s: int) -> np.ndarray:
    """Periodic wrap padding so both dims are divisible by s. Binding rule."""
    x = np.asarray(a, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError(f"Expected 2D array, got shape {x.shape}")

    h, w = x.shape
    pad_h = (-h) % int(s)
    pad_w = (-w) % int(s)

    if pad_h == 0 and pad_w == 0:
        return np.asarray(x, dtype=np.float64, order="C")

    x_pad = np.pad(
        x,
        pad_width=((0, pad_h), (0, pad_w)),
        mode="wrap",
    )
    return np.asarray(x_pad, dtype=np.float64, order="C")


def coarse_grain_block_mean(A: np.ndarray, s: int) -> np.ndarray:
    """Binding coarse graining for one scale s."""
    s = int(s)
    if s <= 0:
        raise ValueError("Scale s must be positive")

    x = np.asarray(A, dtype=np.float64)
    if x.shape != (N_PATCH, N_PATCH):
        raise ValueError(f"Expected shape {(N_PATCH, N_PATCH)}, got {x.shape}")

    x = _wrap_pad_to_multiple(x, s)

    H = x.shape[0] // s
    W = x.shape[1] // s

    # Binding reshape and mean axes (1,3).
    A_s = x.reshape(H, s, W, s).mean(axis=(1, 3))
    return np.asarray(A_s, dtype=np.float64, order="C")


def quantize_u8_patch_internal(x: np.ndarray) -> np.ndarray:
    """Binding quantization to uint8 using per-patch min/max scaling."""
    f = np.asarray(x, dtype=np.float64, order="C")
    xmin = np.min(f)
    xmax = np.max(f)

    if xmax == xmin:
        return np.zeros_like(f, dtype=np.uint8, order="C")

    q_raw = np.floor(np.float64(255.0) * (f - xmin) / (xmax - xmin))
    q_raw = np.clip(q_raw, np.float64(0.0), np.float64(255.0))
    q = q_raw.astype(np.uint8, copy=False)

    return np.asarray(q, dtype=np.uint8, order="C")


def quantized_bytes(q_u8: np.ndarray) -> bytes:
    """Binding byte construction: C-order ravel."""
    q = np.asarray(q_u8, dtype=np.uint8, order="C")
    return q.ravel(order="C").tobytes()
