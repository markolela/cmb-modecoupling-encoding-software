"""CLI entry point for CMB Mode Coupling runs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="cmbmc", add_help=True)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_primary = sub.add_parser("primary", help="Run Sections 2-4 (data + Gaussian null).")
    p_primary.add_argument("--out", required=True, help="Output directory, must be new.")

    p_inj = sub.add_parser("inj", help="Run Section 6 injections.")
    p_inj.add_argument("--out", required=True, help="Output directory, must be new.")
    p_inj.add_argument("--primary-artifacts", required=True, help="Primary artifacts directory (paper/artifacts/primary/<RUN_ID>).")

    p_hm = sub.add_parser("hmdiff", help="Run Section 7 negative control (HMDIFF).")
    p_hm.add_argument("--out", required=True, help="Output directory, must be new.")
    p_hm.add_argument(
        "--primary-artifacts",
        default="paper/artifacts/primary",
        help="Primary artifact directory or its parent (latest subdir is used).",
    )

    p_ffp = sub.add_parser("ffp", help="Run Section 5 interpretative FFP attribution null.")
    p_ffp.add_argument("--out", required=True, help="Output directory, must be new.")
    p_ffp.add_argument("--manifest", default="preregistration/external_inputs/ffp10_manifest_v2.json")
    p_ffp.add_argument("--manifest-ledger", default="preregistration/external_inputs/sha256_ffp10_manifest_v2.txt")
    p_ffp.add_argument(
        "--primary-artifacts",
        default="paper/artifacts/primary",
        help="Primary artifact directory or its parent (latest subdir is used).",
    )

    return ap


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    ap = build_argparser()
    args = ap.parse_args(argv)

    repo = _repo_root()
    out = Path(args.out).expanduser().resolve()
    print("REPO_ROOT:", repo.as_posix())
    print("CMD:", args.cmd)
    print("OUT:", out.as_posix())

    if args.cmd == "primary":
        from cmbmc.run_primary import run_primary
        run_primary(out)
        print("STATUS: primary run completed")
        return 0

    if args.cmd == "ffp":
        man = (repo / str(args.manifest)).resolve()
        led = (repo / str(args.manifest_ledger)).resolve()
        prim = (repo / str(args.primary_artifacts)).resolve()
        print("MANIFEST:", man.as_posix())
        print("MANIFEST_LEDGER:", led.as_posix())
        print("PRIMARY_ARTIFACTS:", prim.as_posix())
        from cmbmc.run_ffp import run_ffp
        run_ffp(out, manifest_path=man, manifest_ledger_path=led, primary_artifacts_dir=prim)
        print("STATUS: ffp run completed")
        return 0

    if args.cmd == "hmdiff":
        prim = (repo / str(args.primary_artifacts)).resolve()
        print("PRIMARY_ARTIFACTS:", prim.as_posix())
        from cmbmc.run_hmdiff import run_hmdiff
        run_hmdiff(out, primary_artifacts_dir=prim)
        print("STATUS: hmdiff run completed")
        return 0
    if args.cmd == "inj":
        prim = (repo / str(args.primary_artifacts)).resolve()
        print("PRIMARY_ARTIFACTS:", prim.as_posix())
        from cmbmc.run_inj import run_inj
        run_inj(out, primary_artifacts_dir=prim)
        print("STATUS: inj run completed")
        return 0


    print("STATUS: not implemented")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
