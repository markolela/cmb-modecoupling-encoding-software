#!/usr/bin/env python3
"""Update paper HMDIFF negative control section from an archived artifact bundle.

Reads:
  paper/artifacts/hmdiff/<RUN_ID>/{hmdiff_summary.json,hmdiff_runinfo.json,hmdiff_metrics.csv,SHA256.txt}

Writes:
  paper/sections/generated/02_results_hmdiff.tex
  paper/tables/hmdiff_effects_table.tex

Validates SHA256.txt against actual files to prevent drift.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REQUIRED = [
    "hmdiff_metrics.csv",
    "hmdiff_runinfo.json",
    "hmdiff_summary.json",
]


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _validate_ledger(bundle: Path) -> None:
    ledger = bundle / "SHA256.txt"
    if not ledger.exists():
        raise SystemExit(f"Missing ledger: {ledger}")

    want = {}
    for line in ledger.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            raise SystemExit(f"Bad ledger line: {line!r}")
        want[parts[1]] = parts[0]

    for name in REQUIRED:
        p = bundle / name
        if not p.exists():
            raise SystemExit(f"Missing artifact: {p}")
        got = _sha256_file(p)
        exp = want.get(name)
        if exp is None:
            raise SystemExit(f"Ledger missing entry for {name}")
        if got != exp:
            raise SystemExit(f"SHA256 mismatch for {name}. GOT={got} EXP={exp}")

    print("OK: SHA256 ledger validated for HMDIFF bundle")


def _load_products(summary_obj: dict) -> dict:
    # Expect SMICA and NILC at top level.
    out = {}
    for prod in ("SMICA", "NILC"):
        if prod in summary_obj and isinstance(summary_obj[prod], dict):
            out[prod] = summary_obj[prod]
    if not out:
        keys = list(summary_obj.keys())
        raise SystemExit(f"Unexpected hmdiff_summary.json format. Keys={keys}")
    return out


def _gate(products: dict) -> bool:
    def lo(prod: str) -> float:
        ci = products[prod]["ci95_delta_theta"]
        return float(ci[0])
    return (lo("SMICA") > 0.0) and (lo("NILC") > 0.0)


def _num(x: float, nd: int = 9) -> str:
    return f"{float(x):.{nd}f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    run_id = str(args.run_id).strip()
    bundle = Path("paper/artifacts/hmdiff") / run_id

    if not bundle.exists():
        raise SystemExit(f"Bundle not found: {bundle}")

    _validate_ledger(bundle)

    summary = json.loads((bundle / "hmdiff_summary.json").read_text(encoding="utf-8"))
    products = _load_products(summary)

    gate = _gate(products)

    # Table.
    tpath = Path("paper/tables/hmdiff_effects_table.tex")
    t = []
    t.append(r"\begin{table}[t]" + "\n")
    t.append(r"\centering" + "\n")
    t.append(r"\caption{HMDIFF negative control (Section 7). The preregistered rule forbids declaring discovery if HMDIFF satisfies the discovery criterion.}" + "\n")
    t.append(r"\label{tab:hmdiff_effects}" + "\n")
    t.append(r"\begin{tabular}{lrrrr}" + "\n")
    t.append(r"\toprule" + "\n")
    t.append(r"Product & $\Delta\theta_{\mathrm{HMDIFF}}$ & CI95 low & CI95 high & $\Delta PP_{\mathrm{HMDIFF}}$ \\" + "\n")
    t.append(r"\midrule" + "\n")
    for prod in ("SMICA", "NILC"):
        x = products[prod]
        lo, hi = x["ci95_delta_theta"]
        dpp = x.get("delta_pp", float("nan"))
        t.append(
            f"{prod} & \\num{{{_num(x['delta_theta'])}}} & "
            f"\\num{{{_num(lo)}}} & \\num{{{_num(hi)}}} & \\num{{{_num(dpp, nd=6)}}} \\\\\n"
        )
    t.append(r"\bottomrule" + "\n")
    t.append(r"\end{tabular}" + "\n")
    t.append(r"\end{table}" + "\n")
    tpath.write_text("".join(t), encoding="utf-8")
    print("WROTE:", tpath.as_posix())

    # Generated section.
    rpath = Path("paper/sections/generated/02_results_hmdiff.tex")
    txt = []
    txt.append("We report the preregistered HMDIFF negative control (Section 7). ")
    txt.append("By protocol, discovery cannot be declared if HMDIFF satisfies the discovery criterion. ")
    if gate:
        txt.append("For this run, HMDIFF yields positive $\\Delta\\theta$ with CI95 fully above zero for both SMICA and NILC. ")
        txt.append("Therefore the preregistered discovery claim is not admissible for this run. \n\n")
    else:
        txt.append("For this run, HMDIFF does not satisfy the discovery criterion. \n\n")

    txt.append(r"\input{tables/hmdiff_effects_table}" + "\n\n")
    txt.append(r"\paragraph{Artifact traceability.}" + "\n")
    txt.append("All numbers in Table~\\ref{tab:hmdiff_effects} are read from ")
    txt.append(r"\path{" + f"paper/artifacts/hmdiff/{run_id}/hmdiff_summary.json" + r"}." + "\n")
    txt.append("File digests for the full HMDIFF artifact bundle are stored in ")
    txt.append(r"\path{" + f"paper/artifacts/hmdiff/{run_id}/SHA256.txt" + r"}." + "\n")
    rpath.write_text("".join(txt), encoding="utf-8")
    print("WROTE:", rpath.as_posix())

    # Ensure 02_results.tex includes the subsection.
    wrapper = Path("paper/sections/02_results.tex")
    w = wrapper.read_text(encoding="utf-8") if wrapper.exists() else ""
    needle = r"\input{sections/generated/02_results_hmdiff}"
    if needle not in w:
        w = w.rstrip() + "\n\n" + r"\subsection{HMDIFF negative control}" + "\n" + needle + "\n"
        wrapper.write_text(w, encoding="utf-8")
        print("UPDATED:", wrapper.as_posix())
    else:
        print("OK: 02_results.tex already references 02_results_hmdiff")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
