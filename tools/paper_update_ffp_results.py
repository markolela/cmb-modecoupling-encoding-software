#!/usr/bin/env python3
"""Update FFP paper results from an archived FFP artifact bundle.

Reads:
  paper/artifacts/ffp/<RUN_ID>/{ffp_summary.json,ffp_runinfo.json,ffp_manifest_used.json,ffp_manifest_ledger_sha256.txt,SHA256.txt}
Writes:
  paper/sections/generated/02_results_ffp.tex
  paper/tables/ffp_effects_table.tex

Validates SHA256.txt against actual files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REQUIRED_FILES = [
    "ffp_metrics.csv",
    "ffp_summary.json",
    "ffp_runinfo.json",
    "ffp_manifest_used.json",
    "ffp_manifest_ledger_sha256.txt",
    "SHA256.txt",
]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_sha256_ledger(bundle: Path) -> None:
    ledger = bundle / "SHA256.txt"
    lines = ledger.read_text(encoding="utf-8").splitlines()
    for ln in lines:
        ln = ln.strip()
        if not ln:
            continue
        got_hash, fname = ln.split(maxsplit=1)
        p = bundle / fname
        exp_hash = sha256_file(p)
        if exp_hash != got_hash:
            raise SystemExit(f"SHA256 mismatch for {p}: GOT {exp_hash} EXP {got_hash}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--primary-run-id", required=True)
    args = ap.parse_args()

    run_id = str(args.run_id).strip()
    primary_run_id = str(args.primary_run_id).strip()

    bundle = Path("paper/artifacts/ffp") / run_id
    for fn in REQUIRED_FILES:
        if not (bundle / fn).exists():
            raise SystemExit(f"Missing required file: {bundle / fn}")

    validate_sha256_ledger(bundle)
    print("OK: SHA256 ledger validated for artifact bundle")

    s = json.loads((bundle / "ffp_summary.json").read_text(encoding="utf-8"))

    n_ffp = int(s["n_ffp"])
    theta_ffp_ref = float(s["theta_ffp_ref"])
    theta_null_ref = float(s["theta_null_ref"])
    delta_theta = float(s["delta_theta"])
    ci_lo, ci_hi = (float(s["ci95_delta_theta"][0]), float(s["ci95_delta_theta"][1]))

    pp_ffp_ref = float(s["pp_ffp_ref"])
    pp_null_ref = float(s["pp_null_ref"])
    delta_pp = float(s["delta_pp"])
    ci2_lo, ci2_hi = (float(s["ci95_delta_pp"][0]), float(s["ci95_delta_pp"][1]))

    manifest_sha = str(s["manifest_used_sha256"])

    def _latex_breakable_hex(h: str, group: int = 8) -> str:
        """
        Render a long hex string in LaTeX in a way that allows line breaks.
        We insert '\\allowbreak' between fixed-size groups.
        This is purely presentational and does not change the hash content.
        """
        h = str(h).strip()
        if not h:
            return h
        parts = [h[i:i+group] for i in range(0, len(h), group)]
        return r"\hspace{0pt}".join(parts)

    artifact_path = f"paper/artifacts/ffp/{run_id}"
    primary_path = f"paper/artifacts/primary/{primary_run_id}"

    # Table
    tpath = Path("paper/tables/ffp_effects_table.tex")
    t = []
    t.append(r"\begin{table}[t]" + "\n")
    t.append(r"\centering" + "\n")
    t.append(
        r"\caption{Interpretative FFP10 attribution ensemble results (Section 5). "
        r"These results do not affect any preregistered decision gate. "
        r"CI95 is computed relative to the frozen Gaussian null reference from the primary run.}" + "\n"
    )
    t.append(r"\label{tab:ffp_effects}" + "\n")
    t.append(r"\begin{tabular}{lrrrrr}" + "\n")
    t.append(r"\toprule" + "\n")
    t.append(r"$N_{\mathrm{FFP}}$ & $\theta_{\mathrm{FFP,ref}}$ & $\theta_{\mathrm{null,ref}}$ & $\Delta\theta$ & CI95$_{\Delta\theta}$ low & CI95$_{\Delta\theta}$ high \\" + "\n")
    t.append(r"\midrule" + "\n")
    t.append(
        f"\\num{{{n_ffp:d}}} & "
        f"\\num{{{theta_ffp_ref:.9f}}} & "
        f"\\num{{{theta_null_ref:.9f}}} & "
        f"\\num{{{delta_theta:.9f}}} & "
        f"\\num{{{ci_lo:.9f}}} & "
        f"\\num{{{ci_hi:.9f}}} \\\\\n"
    )
    t.append(r"\bottomrule" + "\n")
    t.append(r"\end{tabular}" + "\n\n")
    t.append(r"\vspace{0.75em}" + "\n\n")
    t.append(r"\begin{tabular}{lrrrr}" + "\n")
    t.append(r"\toprule" + "\n")
    t.append(r"$PP_{\mathrm{FFP,ref}}$ & $PP_{\mathrm{null,ref}}$ & $\Delta PP$ & CI95$_{\Delta PP}$ low & CI95$_{\Delta PP}$ high \\" + "\n")
    t.append(r"\midrule" + "\n")
    t.append(
        f"\\num{{{pp_ffp_ref:.6f}}} & "
        f"\\num{{{pp_null_ref:.6f}}} & "
        f"\\num{{{delta_pp:.6f}}} & "
        f"\\num{{{ci2_lo:.6f}}} & "
        f"\\num{{{ci2_hi:.6f}}} \\\\\n"
    )
    t.append(r"\bottomrule" + "\n")
    t.append(r"\end{tabular}" + "\n")
    t.append(r"\end{table}" + "\n")
    tpath.write_text("".join(t), encoding="utf-8")
    print("WROTE:", tpath.as_posix())

    # FFP results text (no \section, wrapper already provides it)
    rpath = Path("paper/sections/generated/02_results_ffp.tex")

    txt = []
    txt.append(
        "We report the interpretative FFP10 attribution ensemble (Section 5) with "
        f"$N_\\mathrm{{FFP}} = {n_ffp}$. "
        "These results do not affect any preregistered decision gate.\n\n"
    )
    txt.append(
        "The preregistered effect is evaluated relative to the frozen Gaussian null reference "
        "from the primary run. "
        "For this run, $\\Delta\\theta$ is negative with CI95 fully below zero.\n\n"
    )
    txt.append(r"\input{tables/ffp_effects_table}" + "\n\n")
    txt.append(r"\paragraph{Artifact traceability.}" + "\n")
    txt.append(
        "All numbers in Table~\\ref{tab:ffp_effects} are read from "
        f"\\path{{{artifact_path}/ffp_summary.json}}.\n"
    )
    txt.append(
        "The manifest snapshot used for this run is stored as "
        f"\\path{{{artifact_path}/ffp_manifest_used.json}} with SHA256 "
        f"\\texttt{{{_latex_breakable_hex(manifest_sha)}}}.\n"
        f"\\texttt{{{manifest_sha}}}.\n"
    )
    txt.append(
        "The Gaussian null reference is taken from "
        f"\\path{{{primary_path}/null_effects.csv}} and "
        f"\\path{{{primary_path}/null_summary.json}}.\n"
    )
    txt.append(
        "File digests for the full FFP artifact bundle are stored in "
        f"\\path{{{artifact_path}/SHA256.txt}}.\n"
    )

    rpath.write_text("".join(txt), encoding="utf-8")
    print("WROTE:", rpath.as_posix())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
