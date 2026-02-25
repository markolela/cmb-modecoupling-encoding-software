# scripts/harmonize_beam_nside.py
# Harmonisiert Planck-Karten auf gemeinsamen Beam (FWHM, arcmin) und NSIDE.
# Aufruf:
#   python scripts\harmonize_beam_nside.py --in fits1.fits fits2.fits --fwhm-arcmin 10 --nside 2048 --outdir data\raw\astro\planck\harmonized
from pathlib import Path
import argparse, numpy as np
import healpy as hp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inputs", nargs="+", required=True)
    ap.add_argument("--fwhm-arcmin", type=float, default=10.0)
    ap.add_argument("--nside", type=int, default=2048)
    ap.add_argument("--outdir", default=str(Path("data")/"raw"/"astro"/"planck"/"harmonized"))
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    fwhm_rad = np.deg2rad(args.fwhm_arcmin/60.0)

    for src in args.inputs:
        src = Path(src)
        print(f"[i] Lade {src.name}")
        m = hp.read_map(src, field=None, verbose=False)  # I (oder IQU -> field=None lädt alle)
        # Glätten
        if isinstance(m, np.ndarray) and m.ndim==2:  # IQU
            m_s = []
            for i in range(m.shape[0]):
                m_s.append(hp.sphtfunc.smoothing(m[i], fwhm=fwhm_rad, verbose=False))
            m = np.asarray(m_s)
        else:
            m = hp.sphtfunc.smoothing(m, fwhm=fwhm_rad, verbose=False)
        # NSIDE anpassen
        if isinstance(m, np.ndarray) and m.ndim==2:
            m2 = []
            for i in range(m.shape[0]):
                m2.append(hp.ud_grade(m[i], nside_out=args.nside, pess=False, order_in="RING", order_out="RING"))
            m = np.asarray(m2)
        else:
            m = hp.ud_grade(m, nside_out=args.nside, pess=False, order_in="RING", order_out="RING")

        out = outdir / f"{src.stem}_beam{int(args.fwhm_arcmin)}arcmin_nside{args.nside}.fits"
        hp.write_map(out, m, overwrite=True, dtype=np.float32)
        print(f"[ok] Gespeichert: {out}")

if __name__ == "__main__":
    main()
