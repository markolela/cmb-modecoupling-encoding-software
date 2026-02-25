# scripts/smooth_patches_planar.py
# Angleicht Patches eines Datensatzes auf eine gemeinsame Ziel-FWHM (2D-Gauss) – ohne healpy.
# Annahme: kleine FoVs (z.B. 12°) ⇒ planare Gauss-Glättung ist hinreichend.

import json, math, shutil
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
PATCHES = REPO / "data" / "processed" / "astro" / "patches"

def fft_gauss_blur(img, sigma_px: float):
    if sigma_px <= 0: 
        return img
    H, W = img.shape
    fy = np.fft.fftfreq(H)
    fx = np.fft.fftfreq(W)
    FY, FX = np.meshgrid(fy, fx, indexing="ij")
    G = np.exp(-2.0 * (np.pi**2) * (sigma_px**2) * (FX**2 + FY**2))
    F = np.fft.fft2(img)
    return np.fft.ifft2(F * G).real

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="dataset_tag (Ordner unter patches/)")
    ap.add_argument("--target-fwhm-arcmin", type=float, required=True)
    ap.add_argument("--assume-native-fwhm-arcmin", type=float, default=5.0,
                    help="Annahme über native Beam-FWHM in arcmin (Planck CMB ~5')")
    ap.add_argument("--out-suffix", default=None, help="z.B. _sm10am (für 10 arcmin)")
    args = ap.parse_args()

    src = PATCHES / args.dataset
    if not src.exists():
        raise SystemExit(f"Dataset-Ordner fehlt: {src}")

    meta_path = src / "meta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    N = int(meta.get("N", 1024))
    fov_deg = float(meta.get("fov_deg", 12.0))

    # Pixelmaßstab
    pix_arcmin = (fov_deg * 60.0) / N
    # Sigma in Pixel
    fwhm2sigma = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    sigma_tgt_px = (args.target_fwhm_arcmin / pix_arcmin) * fwhm2sigma
    sigma_nat_px = (args.assume_native_fwhm_arcmin / pix_arcmin) * fwhm2sigma
    delta_sigma_px = max(0.0, math.sqrt(max(0.0, sigma_tgt_px**2 - sigma_nat_px**2)))

    out_suffix = args.out_suffix or f"_sm{int(round(args.target_fwhm_arcmin))}am"
    dst_tag = args.dataset + out_suffix
    dst = PATCHES / dst_tag
    dst.mkdir(parents=True, exist_ok=True)

    # Patches glätten
    patch_files = sorted(src.glob("patch_*.npy"))
    smoothed = []
    for pf in patch_files:
        arr = np.load(pf)
        arr_sm = fft_gauss_blur(arr, delta_sigma_px)
        np.save(dst / pf.name, arr_sm.astype(arr.dtype))
        smoothed.append(arr_sm)

    # Stack neu schreiben
    stack_path = dst / f"{dst_tag}_stack.npz"
    if smoothed:
        np.savez_compressed(stack_path, patches=np.stack(smoothed, axis=0))

    # Meta copy and update so that meta.json points to the smoothed patches
    meta_out = dict(meta)
    meta_out["dataset"] = dst_tag
    meta_out["derived_from"] = args.dataset
    meta_out["patches"] = [str(dst / pf.name) for pf in patch_files]
    meta_out["stack"] = str(stack_path)

    meta_out["smoothing"] = {
        "type": "gaussian_planar",
        "target_fwhm_arcmin": args.target_fwhm_arcmin,
        "assumed_native_fwhm_arcmin": args.assume_native_fwhm_arcmin,
        "delta_sigma_px": delta_sigma_px,
        "pix_arcmin": pix_arcmin,
    }
    (dst / "meta.json").write_text(json.dumps(meta_out, indent=2), encoding="utf-8")

    print(f"[OK] Glättung fertig: {dst_tag}")
    print(f"    Ziel-FWHM={args.target_fwhm_arcmin}′, angenommen native={args.assume_native_fwhm_arcmin}′")
    print(f"    Δσ(px)={delta_sigma_px:.3f}, Pixelmaßstab={pix_arcmin:.3f}′/px")

if __name__ == "__main__":
    main()
