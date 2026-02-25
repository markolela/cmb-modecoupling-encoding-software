# scripts/jackknife_hemi_t3.py  (fix: robustes Auslesen von b aus meta["centers"])
import json, subprocess, sys
from pathlib import Path
import argparse
import numpy as np
import re
import shutil

REPO   = Path(__file__).resolve().parents[1]
PATCHD = REPO / "data" / "processed" / "astro" / "patches"

def read_meta(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))

def _parse_min_abs_b_from_tag(tag: str):
    m = re.search(r"_b(\d+)_", tag)
    return float(m.group(1)) if m else None

def extract_b_list(meta: dict, tag: str):
    """
    Robust: bevorzugt meta['patches']-Einträge mit b/lat, ansonsten meta['centers'].
    Danach auf die tatsächlich verwendete Teilmenge abbilden:
      - wenn meta['subset']['indices'] existiert -> damit indizieren
      - sonst Fallback: über |b|>=bmin aus dem Tag (z.B. _b65_), falls vorhanden
    """
    # ---- (A) Alle b-Werte sammeln (wie bisher) ----
    pl = meta.get("patches")
    if isinstance(pl, list) and pl and isinstance(pl[0], dict):
        b_all = []
        for p in pl:
            for key in ("b_deg","lat_deg","gal_b_deg","galactic_b_deg","glat_deg"):
                if key in p:
                    b_all.append(float(p[key])); break
            else:
                if "center_gal_deg" in p and isinstance(p["center_gal_deg"], (list, tuple)) and len(p["center_gal_deg"])>=2:
                    b_all.append(float(p["center_gal_deg"][1]))
                elif "center" in p and isinstance(p["center"], (list, tuple)) and len(p["center"])>=2 and str(p.get("coord","galactic")).lower().startswith("gal"):
                    b_all.append(float(p["center"][1]))
                else:
                    raise KeyError("Patch-Eintrag ohne b/lat-Feld gefunden.")
    else:
        centers = meta.get("centers")
        if isinstance(centers, list) and centers and isinstance(centers[0], dict):
            frame = str(meta.get("frame","")).lower()
            if "gal" in frame and "lat_deg" in centers[0]:
                b_all = [float(c["lat_deg"]) for c in centers]
            else:
                raise KeyError("centers vorhanden, aber frame nicht galaktisch oder lat_deg fehlt.")
        else:
            for k in ("centers_gal_deg","centers_galactic_deg","centers_gal_lonlat_deg","centers_deg","centers_lonlat_deg"):
                v = meta.get(k)
                if isinstance(v, list) and v and isinstance(v[0], (list,tuple)) and len(v[0])>=2:
                    b_all = [float(p[1]) for p in v]
                    break
            else:
                raise KeyError("Konnte galaktische Breiten (b) nicht aus meta.json extrahieren.")

    # ---- (B) Auf verwendete Teilmenge einschränken ----
    sub = meta.get("subset") or {}
    idx = sub.get("indices")
    if isinstance(idx, list) and idx and isinstance(idx[0], int):
        return [b_all[i] for i in idx]

    # Fallback: |b|-Schwelle aus dem Datensatz-Tag
    bmin = _parse_min_abs_b_from_tag(tag)
    if bmin is not None:
        return [b for b in b_all if abs(b) >= bmin]

    # Wenn keine Teilmenge ermittelbar: alles zurückgeben (kann zu Mismatch führen)
    return b_all

def subset_and_write(src_tag: str, hemi: str):
    src = PATCHD / src_tag
    if not src.exists():
        raise FileNotFoundError(f"Patch-Ordner fehlt: {src}")
    meta = read_meta(src / "meta.json")
    b_list = extract_b_list(meta, src_tag)
    patch_files = sorted(
        src.glob("patch_*.npy"),
        key=lambda p: int(re.search(r"patch_(\d+)", p.stem).group(1)) if re.search(r"patch_(\d+)", p.stem) else p.name
    )

    # Debug-Hinweis
    if len(patch_files) != len(b_list):
        print(f"[DBG] src_tag={src_tag}  subset={meta.get('subset',{})}  parsed_bmin={_parse_min_abs_b_from_tag(src_tag)}")

    if len(patch_files) != len(b_list):
        raise RuntimeError(f"Anzahl Patches ({len(patch_files)}) != Anzahl b-Werte ({len(b_list)})")

    keep_idx = [i for i,b in enumerate(b_list) if (b>=0.0) == (hemi=="N")]
    if not keep_idx:
        raise RuntimeError(f"Hemisphäre {hemi}: keine Patches ausgewählt.")

    dst_tag = f"{src_tag}_hemi{hemi}"
    dst = PATCHD / dst_tag
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)

    stack=[]
    for i in keep_idx:
        arr = np.load(patch_files[i]); stack.append(arr)
        np.save(dst / patch_files[i].name, arr)
    np.savez_compressed(dst / f"{dst_tag}_stack.npz", patches=np.stack(stack, axis=0))

    meta_out = dict(meta)
    meta_out["dataset"] = dst_tag
    meta_out["parent_tag"] = src_tag
    meta_out["subset"] = {
        "type": "hemisphere",
        "rule": "gal_b>=0" if hemi == "N" else "gal_b<0",
        "parent_indices": keep_idx,
        "count": len(keep_idx),
        "source": "jackknife_hemi_t3",
    }
    meta_out["hemi"] = hemi
    meta_out["mask_label"] = meta.get("mask_label", f"|b|≥{int(_parse_min_abs_b_from_tag(src_tag) or 0)}; FOV 12°")
    meta_out["created_at"] = __import__("datetime").datetime.now(
        __import__("datetime").timezone.utc
    ).isoformat(timespec="seconds")
    meta_out["subset"]["source"] = "jackknife_hemi_t3"

    (dst / "meta.json").write_text(json.dumps(meta_out, indent=2), encoding="utf-8")
    return dst_tag

def run_t3(tag: str, trend: str, scales: str, null_n: int, null_seed: int):
    cmd = [
        sys.executable, "-m", "scripts.run_t3_on_patches",
        "--dataset", tag,
        "--scales", scales,
        "--trend", trend,
        "--agg", "median",
        "--jobs-data", "20",
        "--null", str(null_n),
        "--null-seed", str(null_seed),
        "--null-family", "phase_randomized",
        "--jobs-null", "20",
    ]
    print(">>", " ".join(cmd))
    subprocess.check_call(cmd, cwd=REPO)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", required=True, help="Kommagetrennt: <tag1>,<tag2>,...")
    ap.add_argument("--trend", default="log")
    ap.add_argument("--scales", default="1,2,4")
    ap.add_argument("--null-n", type=int, default=200)
    ap.add_argument("--null-seed", type=int, default=12345)
    args = ap.parse_args()

    tags = [t.strip() for t in args.datasets.split(",") if t.strip()]
    for base_tag in tags:
        print(f"[INFO] Jackknife für: {base_tag}")
        n_tag = subset_and_write(base_tag, "N")
        s_tag = subset_and_write(base_tag, "S")
        for t in (n_tag, s_tag):
            run_t3(t, args.trend, args.scales, args.null_n, args.null_seed)

if __name__ == "__main__":
    main()
