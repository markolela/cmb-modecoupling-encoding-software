# scripts/harmonize_planck_hm_sm10am.py
# Bringt Planck HM1/HM2 CMB I-Map auf 10′ FWHM und schreibt nach .../planck/harmonized/.
# Benennung exakt wie im Patch-Skript erwartet.

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np

import healpy as hp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="src",  required=True, help="Input FITS (Planck COM_CMB_IQU-..._halfmission-*.fits)")
    ap.add_argument("--out", dest="dest", required=True, help="Output FITS (…_sm10am.fits in .../planck/harmonized/)")
    args = ap.parse_args()

    src  = Path(args.src)
    dest = Path(args.dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Feld 0 = I/T
    print(f"[READ] {src}")
    m = hp.read_map(src.as_posix(), field=0, verbose=False)
    m = m.astype(np.float32, copy=False)

    # 10 arcmin FWHM in Radiant
    fwhm = np.deg2rad(10.0/60.0)
    print(f"[SMOOTH] FWHM = 10′")
    ms = hp.smoothing(m, fwhm=fwhm, verbose=False).astype(np.float32, copy=False)

    print(f"[WRITE] {dest}")
    hp.write_map(dest.as_posix(), ms, dtype=np.float32, overwrite=True)

    print("[OK] Done.")

if __name__ == "__main__":
    main()
