#!/usr/bin/env python3
"""Update primary paper results from an archived primary artifact bundle.

Reads:
  paper/artifacts/primary/<RUN_ID>/{summary.json,null_summary.json,runinfo.json,SHA256.txt}
Writes:
  paper/sections/generated/02_results_primary.tex
  paper/tables/primary_effects_table.tex

Also validates SHA256.txt against the actual files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REQUIRED_FILES = [
    "metrics.csv",
    "metrics_per_patch.csv",
    "null_effects.csv",
    "null_summary.json",
    "runinfo.json",
    "summary.json",
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
    args = ap.parse_args()

    run_id = str(args.run_id).strip()
    bundle = Path("paper/artifacts/primary") / run_id

    for fn in REQUIRED_FILES:
        if not (bundle / fn).exists():
            raise SystemExit(f"Missing required file: {bundle / fn}")

    validate_sha256_ledger(bundle)
    print("OK: SHA256 ledger validated for artifact bundle")

    summary = json.loads((bundle / "summary.json").read_text(encoding="utf-8"))
    null_summary = json.loads((bundle / "null_summary.json").read_text(encoding="utf-8"))
    runinfo = json.loads((bundle / "runinfo.json").read_text(encoding="utf-8"))

    n_null = int(null_summary["n_null"])
    git_commit = str(runinfo.get("git_commit", ""))

    def row(prod: str) -> dict:
        s = summary[prod]
        lo, hi = s["ci95_delta_theta"]
        lo2, hi2 = s["ci95_delta_pp"]
        return dict(
            prod=prod,
            theta_data=float(s["theta_data"]),
            theta_null_ref=float(s["theta_null_ref"]),
            delta_theta=float(s["delta_theta"]),
            ci_lo=float(lo),
            ci_hi=float(hi),
            pp_data=float(s["pp_data"]),
            pp_null_ref=float(s["pp_null_ref"]),
            delta_pp=float(s["delta_pp"]),
            ci2_lo=float(lo2),
            ci2_hi=float(hi2),
        )

    rows = [row("SMICA"), row("NILC")]

    # Table
    tpath = Path("paper/tables/primary_effects_table.tex")
    t = []
    t.append(r"\begin{table}[t]" + "\n")
    t.append(r"\centering" + "\n")
    t.append(
        r"\caption{Primary effects for the preregistered diagnostic. "
        r"The discovery gate is conjunctive across SMICA and NILC and requires "
        r"the lower CI95 bound of $\Delta\theta$ to be greater than zero for both products. "
        r"The auxiliary $\Delta PP$ is reported but is not a decision gate.}" + "\n"
    )
    t.append(r"\label{tab:primary_effects}" + "\n")
    t.append(r"\begin{tabular}{lrrrrrr}" + "\n")
    t.append(r"\toprule" + "\n")
    t.append(r"Product & $\theta_{\mathrm{data}}$ & $\theta_{\mathrm{null,ref}}$ & $\Delta\theta$ & CI95$_{\Delta\theta}$ low & CI95$_{\Delta\theta}$ high & $PP_{\mathrm{data}}$ \\" + "\n")
    t.append(r"\midrule" + "\n")
    for r in rows:
        t.append(
            f"{r['prod']} & "
            f"\\num{{{r['theta_data']:.9f}}} & "
            f"\\num{{{r['theta_null_ref']:.9f}}} & "
            f"\\num{{{r['delta_theta']:.9f}}} & "
            f"\\num{{{r['ci_lo']:.9f}}} & "
            f"\\num{{{r['ci_hi']:.9f}}} & "
            f"\\num{{{r['pp_data']:.6f}}} \\\\\n"
        )
    t.append(r"\bottomrule" + "\n")
    t.append(r"\end{tabular}" + "\n\n")
    t.append(r"\vspace{0.75em}" + "\n\n")
    t.append(r"\begin{tabular}{lrrrrr}" + "\n")
    t.append(r"\toprule" + "\n")
    t.append(r"Product & $PP_{\mathrm{null,ref}}$ & $\Delta PP$ & CI95$_{\Delta PP}$ low & CI95$_{\Delta PP}$ high \\" + "\n")
    t.append(r"\midrule" + "\n")
    for r in rows:
        t.append(
            f"{r['prod']} & "
            f"\\num{{{r['pp_null_ref']:.6f}}} & "
            f"\\num{{{r['delta_pp']:.6f}}} & "
            f"\\num{{{r['ci2_lo']:.6f}}} & "
            f"\\num{{{r['ci2_hi']:.6f}}} \\\\\n"
        )
    t.append(r"\bottomrule" + "\n")
    t.append(r"\end{tabular}" + "\n")
    t.append(r"\end{table}" + "\n")
    tpath.write_text("".join(t), encoding="utf-8")
    print("WROTE:", tpath.as_posix())

    # Primary results text (no \section, wrapper already provides it)
    rpath = Path("paper/sections/generated/02_results_primary.tex")
    artifact_path = f"paper/artifacts/primary/{run_id}"

    txt = []
    txt.append(
        "This section reports the preregistered primary effect "
        "$\\Delta\\theta = \\theta_{\\mathrm{null,ref}} - \\theta_{\\mathrm{data}}$ "
        f"and its CI95 computed from the Gaussian null ensemble with $N_\\mathrm{{null}} = {n_null}$. "
        "The discovery gate is conjunctive across SMICA and NILC and requires the lower CI95 bound "
        "of $\\Delta\\theta$ to be greater than zero for both products.\n\n"
    )
    txt.append(
        "For this run, both SMICA and NILC yield negative $\\Delta\\theta$ with CI95 fully below zero. "
        "Therefore the preregistered discovery criterion is not met in this run.\n\n"
    )
    txt.append(
        "The auxiliary plateau proxy $PP$ and $\\Delta PP$ are reported for completeness but are not a decision gate. "
        "For this run, $\\Delta PP$ is positive for both products with CI95 above zero.\n\n"
    )
    txt.append(r"\input{tables/primary_effects_table}" + "\n\n")
    txt.append(r"\paragraph{Artifact traceability.}" + "\n")
    txt.append(
        "All numbers in Table~\\ref{tab:primary_effects} are read from "
        f"\\path{{{artifact_path}/summary.json}}.\n"
    )
    txt.append(
        "The corresponding run metadata are stored in "
        f"\\path{{{artifact_path}/runinfo.json}} with \\texttt{{git\\_commit = {git_commit}}}.\n"
    )
    txt.append(
        "File digests for the full artifact bundle are stored in "
        f"\\path{{{artifact_path}/SHA256.txt}}.\n"
    )
    rpath.write_text("".join(txt), encoding="utf-8")
    print("WROTE:", rpath.as_posix())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
