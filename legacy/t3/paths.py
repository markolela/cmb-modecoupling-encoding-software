# scripts/t3/paths.py
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict


# Wichtig.
# Diese Datei liegt unter scripts/t3.
# Damit REPO identisch bleibt wie vorher in scripts/run_t3_on_patches.py,
# muss hier parents[2] verwendet werden.
REPO = Path(__file__).resolve().parents[2]
PATCH = REPO / "data" / "processed" / "astro" / "patches"
OUT = REPO / "data" / "processed" / "astro" / "suite"


def load_manifest() -> Dict[str, dict]:
    mani = PATCH / "patches_manifest.json"
    if mani.exists():
        return json.loads(mani.read_text(encoding="utf-8"))
    return {}


def load_dataset_meta(ds: str) -> dict:
    meta_path = PATCH / ds / "meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing meta file: {meta_path}")
    return json.loads(meta_path.read_text(encoding="utf-8"))



def resolve_patch_paths(ds: str, meta: dict) -> list[Path]:
    """
    Robust patch file path resolution.
    Priority:
    1) Use meta["patches"] paths if they exist.
    2) If absolute paths do not exist, try PATCH/ds/<basename>.
    3) Fallback: PATCH/ds/patch_*.npy sorted.
    """
    ds_dir = PATCH / ds

    raw = list(meta.get("patches", []))
    out: list[Path] = []

    for p in raw:
        pp = Path(p)
        if pp.is_absolute():
            if pp.exists():
                out.append(pp)
                continue
            cand = ds_dir / pp.name
            if cand.exists():
                out.append(cand)
                continue
            out.append(pp)
            continue

        cand = REPO / pp
        if cand.exists():
            out.append(cand)
            continue

        cand2 = ds_dir / pp
        if cand2.exists():
            out.append(cand2)
            continue

        out.append(cand)

    if out and all(p.exists() for p in out):
        return out

    fallback = sorted(ds_dir.glob("patch_*.npy"))
    if not fallback:
        raise FileNotFoundError(f"No patch files found under {ds_dir}")
    return fallback
