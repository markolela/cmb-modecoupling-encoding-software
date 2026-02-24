# -*- coding: utf-8 -*-
"""
Generate a deterministic, versioned, frozen patch-center list for the CMB pipeline.

This script is protocol-facing. It must remain stable once v1 is frozen.
It uses only the fixed mask product to avoid centers dominated by masked zeros.

Output:
- CSV: l_deg, b_deg, theta_rad, phi_rad, cand_pix, unmasked_fraction
- JSON manifest: parameters, mask hash, selection rules
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Tuple

import numpy as np
import healpy as hp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def angdist_rad(v1: np.ndarray, v2: np.ndarray) -> float:
    # v1 and v2 are 3-vectors on the unit sphere
    c = float(np.clip(np.dot(v1, v2), -1.0, 1.0))
    return float(np.arccos(c))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mask", required=True, type=Path, help="Path to fixed mask FITS.")
    ap.add_argument("--out_csv", required=True, type=Path)
    ap.add_argument("--out_manifest", required=True, type=Path)
    ap.add_argument("--n_centers", type=int, default=256)
    ap.add_argument("--nside_cand", type=int, default=64)
    ap.add_argument("--radius_deg", type=float, default=8.5)
    ap.add_argument("--min_unmasked", type=float, default=0.95)
    ap.add_argument("--center_threshold", type=float, default=0.999)
    args = ap.parse_args()

    mask_path: Path = args.mask
    if not mask_path.exists():
        raise FileNotFoundError(mask_path)

    mask_hash = sha256_file(mask_path)

    # Read mask as a HEALPix map
    mask_map = hp.read_map(str(mask_path), verbose=False)
    nside_mask = hp.get_nside(mask_map)

    # Candidate centers from HEALPix pixel centers
    nside_cand = int(args.nside_cand)
    cand_pix = np.arange(hp.nside2npix(nside_cand), dtype=np.int64)
    theta, phi = hp.pix2ang(nside_cand, cand_pix, nest=False)
    cand_vecs = hp.ang2vec(theta, phi)  # shape (N, 3)

    radius_rad = np.deg2rad(float(args.radius_deg))
    min_unmasked = float(args.min_unmasked)
    center_thr = float(args.center_threshold)

    # Compute unmasked fraction in a disc around each candidate
    unmasked_frac = np.zeros(cand_pix.shape[0], dtype=np.float64)
    center_ok = np.zeros(cand_pix.shape[0], dtype=bool)

    for i in range(cand_pix.shape[0]):
        v = cand_vecs[i]
        disc = hp.query_disc(nside_mask, v, radius_rad, inclusive=False, nest=False)
        vals = mask_map[disc]
        # Treat the mask as effectively binary by thresholding near 1
        unmasked_frac[i] = float(np.mean(vals >= center_thr))
        # Center pixel must be unmasked
        t0, p0 = float(theta[i]), float(phi[i])
        pix0 = int(hp.ang2pix(nside_mask, t0, p0, nest=False))
        center_ok[i] = bool(mask_map[pix0] >= center_thr)

    keep = np.where((center_ok) & (unmasked_frac >= min_unmasked))[0]
    if keep.size < int(args.n_centers):
        raise RuntimeError(
            f"Not enough candidates after mask filtering. "
            f"Kept {keep.size}, need {args.n_centers}."
        )

    keep_pix = cand_pix[keep]
    keep_theta = theta[keep]
    keep_phi = phi[keep]
    keep_vecs = cand_vecs[keep]
    keep_frac = unmasked_frac[keep]

    # Greedy farthest-point selection with deterministic tie breaks
    n_centers = int(args.n_centers)

    # Start: max unmasked_fraction, tie break by smallest pixel id
    start_idx = np.lexsort((keep_pix, -keep_frac))[0]
    selected = [int(start_idx)]
    selected_vecs = [keep_vecs[start_idx]]

    # Precompute for speed
    keep_vecs_np = np.asarray(keep_vecs, dtype=np.float64)

    for _ in range(1, n_centers):
        best_i = None
        best_min_d = -1.0

        for i in range(keep_vecs_np.shape[0]):
            if i in selected:
                continue
            v = keep_vecs_np[i]
            # min distance to selected set
            min_d = min(angdist_rad(v, sv) for sv in selected_vecs)

            if min_d > best_min_d:
                best_min_d = min_d
                best_i = i
            elif min_d == best_min_d:
                # tie break by smallest keep_pix
                if best_i is not None and int(keep_pix[i]) < int(keep_pix[best_i]):
                    best_i = i

        if best_i is None:
            raise RuntimeError("Selection failed unexpectedly.")

        selected.append(int(best_i))
        selected_vecs.append(keep_vecs_np[best_i])

    sel_pix = keep_pix[selected]
    sel_theta = keep_theta[selected]
    sel_phi = keep_phi[selected]
    sel_frac = keep_frac[selected]

    # Convert to Galactic lon, lat as degrees, assuming Galactic HEALPix convention
    l_deg = np.rad2deg(sel_phi) % 360.0
    b_deg = 90.0 - np.rad2deg(sel_theta)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)

    # Write CSV
    header = "idx,l_deg,b_deg,theta_rad,phi_rad,cand_pix,unmasked_fraction"
    rows = []
    for k in range(n_centers):
        rows.append(
            f"{k},"
            f"{l_deg[k]:.6f},"
            f"{b_deg[k]:.6f},"
            f"{float(sel_theta[k]):.12f},"
            f"{float(sel_phi[k]):.12f},"
            f"{int(sel_pix[k])},"
            f"{float(sel_frac[k]):.6f}"
        )
    args.out_csv.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")

    # Write manifest
    manifest = {
        "version": "patch_centers_n256_fov12_galactic_v1",
        "coord_frame": "galactic",
        "fov_deg": 12.0,
        "patch_grid": [256, 256],
        "mask_path": str(mask_path),
        "mask_sha256": mask_hash,
        "mask_center_threshold": center_thr,
        "nside_mask": int(nside_mask),
        "nside_candidate": int(nside_cand),
        "candidate_order": "RING, pix ascending",
        "disc_radius_deg": float(args.radius_deg),
        "min_unmasked_fraction": float(min_unmasked),
        "selection": "greedy_farthest_point_on_sphere",
        "tie_break": "smallest_candidate_pix",
        "n_centers": int(n_centers),
        "outputs": {
            "csv": str(args.out_csv),
        },
    }
    args.out_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()