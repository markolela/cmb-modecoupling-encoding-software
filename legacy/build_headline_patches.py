# build_headline_patches.py


from __future__ import annotations

import argparse
from pathlib import Path

import yaml  # pip install pyyaml

# Import the existing patch builder
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from make_real_patches_standalone import build_dataset  # type: ignore

RAW_PLANCK = REPO / "data" / "raw" / "astro" / "planck"
RAW_PLANCK_HARM = RAW_PLANCK / "harmonized"

FITS_BY_MAP = {
    "SMICA": RAW_PLANCK / "COM_CMB_IQU-smica_2048_R3.00_full.fits",
    "NILC":  RAW_PLANCK / "COM_CMB_IQU-nilc_2048_R3.00_full.fits",
    "SMICA_HM1": RAW_PLANCK_HARM / "COM_CMB_IQU-smica_2048_R3.00_hm1_sm10am.fits",
    "SMICA_HM2": RAW_PLANCK_HARM / "COM_CMB_IQU-smica_2048_R3.00_hm2_sm10am.fits",
    "NILC_HM1":  RAW_PLANCK_HARM / "COM_CMB_IQU-nilc_2048_R3.00_hm1_sm10am.fits",
    "NILC_HM2":  RAW_PLANCK_HARM / "COM_CMB_IQU-nilc_2048_R3.00_hm2_sm10am.fits",
}

N_PIX = 1024  # fixed patch resolution


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/headline.yaml")
    args = ap.parse_args()

    cfg_path = (REPO / args.config).resolve()
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))

    fov = float(cfg["patches"]["fov_deg"])
    bcut = int(cfg["mask"]["b_cut_deg"])
    n_patches = int(cfg["patches"]["n_patches"])
    seed = int(cfg["patches"].get("seed", 41))  # default is fixed and explicit

    sel = str(cfg["patches"].get("selection", "fixed")).lower().strip()
    if sel != "fixed":
        raise SystemExit(f"patches.selection muss 'fixed' sein, got: {sel!r}")

    beam = int(cfg["beam"]["fwhm_arcmin"])
    sm_tag = f"sm{beam}am"
    dataset = str(cfg.get("dataset", "planck2018"))

    ref_map = str(cfg["maps"]["reference"])
    cmp_map = str(cfg["maps"]["comparison"])

    def dataset_tag(map_name: str) -> str:
        return f"{dataset}_{map_name}_T_fov{int(fov)}_b{bcut}_n{n_patches}_{sm_tag}"

    def fits_for(map_name: str) -> Path:
        p = FITS_BY_MAP.get(map_name)
        if p is None:
            raise SystemExit(f"Unknown map name in config: {map_name}. Known: {sorted(FITS_BY_MAP)}")
        if not p.exists():
            raise SystemExit(f"Missing FITS for {map_name}: {p}")
        return p

    for map_name in (ref_map, cmp_map):
        ds = dataset_tag(map_name)
        fits = fits_for(map_name)
        print(f"Building patches: {ds}")
        build_dataset(
            ds,
            fits,
            n_patches=n_patches,
            N=N_PIX,
            fov_deg=fov,
            lat_cut_deg=float(bcut),
            seed=seed,
        )

    print("[OK] Base patches built.")


if __name__ == "__main__":
    main()
