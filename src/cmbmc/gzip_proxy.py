"""Gzip compression proxy and baseline correction, binding Section 2.7.

Binding rules implemented.
- gzip.compress(data_bytes, compresslevel=6, mtime=0)
- compressed_bytes = len(gzip_bytes) including gzip header/footer
- bpc(p,s) = compressed_bytes / number_of_cells
- Baseline bpc0(s) computed once per scale using:
  rng = numpy.random.Generator(numpy.random.PCG64(baseline_seed))
  For each scale s in s_list = [1,2,4]:
    Draw exactly three uint8 arrays sequentially from the same generator.
- Encoding proxy:
  kappa(p,s) = (bpc(p,s) - bpc0(s)) * log(256)
- Canonical gzip smoketest hash must match committed repo-root file gzip_smoketest_hash.txt
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip
import hashlib

import numpy as np


S_LIST = [1, 2, 4]
GZIP_LEVEL = 6
LOG_256 = np.log(np.float64(256.0))


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def gzip_compressed_length(data_bytes: bytes) -> int:
    gz = gzip.compress(data_bytes, compresslevel=int(GZIP_LEVEL), mtime=0)
    return int(len(gz))


def _read_smoketest_ref() -> str:
    p = repo_root() / "gzip_smoketest_hash.txt"
    if not p.exists():
        raise FileNotFoundError(p.as_posix())
    s = p.read_text(encoding="utf-8")
    if not s.endswith("\n"):
        raise ValueError("gzip_smoketest_hash.txt must end with a trailing newline")
    h = s.strip()
    if len(h) != 64 or any(c not in "0123456789abcdef" for c in h):
        raise ValueError("gzip_smoketest_hash.txt must contain lowercase sha256 hex digest only")
    return h


def gzip_smoketest_or_raise() -> None:
    """Binding smoketest."""
    smoke = np.concatenate([np.arange(256, dtype=np.uint8)] * 4)
    smoke_bytes = smoke.tobytes()
    gz = gzip.compress(smoke_bytes, compresslevel=int(GZIP_LEVEL), mtime=0)
    got = hashlib.sha256(gz).hexdigest()
    exp = _read_smoketest_ref()
    if got != exp:
        raise ValueError(f"gzip smoketest mismatch got={got} exp={exp}")


def baseline_bpc0_by_scale(*, baseline_seed: int) -> dict[int, np.float64]:
    """Binding baseline construction. Returns {s: bpc0(s)} for s in S_LIST."""
    rng = np.random.Generator(np.random.PCG64(int(baseline_seed)))

    out: dict[int, np.float64] = {}
    for s in S_LIST:
        H = 256 // int(s)
        W = 256 // int(s)
        number_of_cells = int(H * W)

        vals: list[np.float64] = []
        for _k in range(3):
            U = rng.integers(0, 256, size=(H, W), dtype=np.uint8)
            data_bytes = U.ravel(order="C").tobytes()
            comp_len = gzip_compressed_length(data_bytes)
            vals.append(np.float64(comp_len) / np.float64(number_of_cells))

        out[int(s)] = np.mean(vals, dtype=np.float64)

    return out


def bpc_from_quantized_u8(q_u8: np.ndarray) -> np.float64:
    """Compute bpc for one quantized field q (uint8) using binding gzip rule."""
    q = np.asarray(q_u8, dtype=np.uint8, order="C")
    data_bytes = q.ravel(order="C").tobytes()
    comp_len = gzip_compressed_length(data_bytes)
    number_of_cells = int(q.size)
    return np.float64(comp_len) / np.float64(number_of_cells)


def kappa_from_bpc(bpc: np.float64, bpc0: np.float64) -> np.float64:
    return (np.float64(bpc) - np.float64(bpc0)) * np.float64(LOG_256)
