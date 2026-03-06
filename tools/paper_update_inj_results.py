#!/usr/bin/env python3
"""Update paper injection suite section from an archived artifact bundle.

Reads:
  paper/artifacts/inj/<RUN_ID>/{inj_metrics.csv,inj_null_effects.csv,inj_runinfo.json,SHA256.txt}

Writes:
  paper/sections/generated/02_results_inj.tex
  paper/tables/inj_effects_table.tex

Also ensures paper/sections/02_results.tex includes the generated subsection.

All values are computed deterministically from artifacts:
- theta_inj(alpha) and PP_inj(alpha) are medians over injection realizations in inj_metrics.csv
- theta_null and PP_null distributions are read from inj_null_effects.csv
- CI95 uses numpy.quantile(..., method="linear") per plan.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


REQ = ["inj_metrics.csv", "inj_null_effects.csv", "inj_runinfo.json"]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_sha256_ledger(bundle_dir: Path) -> None:
    ledger = bundle_dir / "SHA256.txt"
    want: dict[str, str] = {}
    for line in ledger.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            raise SystemExit(f"Bad SHA256.txt line: {line!r}")
        want[parts[1]] = parts[0]

    for name in REQ:
        if name not in want:
            raise SystemExit(f"SHA256.txt missing entry for {name}")
        got = sha256_file(bundle_dir / name)
        if got != want[name]:
            raise SystemExit(f"SHA256 mismatch for {name}: got {got}, expected {want[name]}")

    print("OK: SHA256 ledger validated for injection bundle")


def _pick_col(fieldnames: list[str], candidates: list[str]) -> str:
    s = set(fieldnames)
    for c in candidates:
        if c in s:
            return c
    raise SystemExit(f"Could not find any of columns {candidates} in header={fieldnames}")


def read_null_arrays(null_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    with null_csv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        if r.fieldnames is None:
            raise SystemExit("null_effects.csv has no header")
        theta_col = _pick_col(r.fieldnames, ["theta", "theta_null", "theta_r", "theta_null_r"])
        pp_col = _pick_col(r.fieldnames, ["pp", "PP", "pp_null", "pp_r", "PP_r", "pp_null_r", "PP_null"])
        theta: list[float] = []
        pp: list[float] = []
        for row in r:
            theta.append(float(row[theta_col]))
            pp.append(float(row[pp_col]))
    return np.asarray(theta, dtype=np.float64), np.asarray(pp, dtype=np.float64)


def read_inj_by_alpha(inj_csv: Path) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    with inj_csv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        if r.fieldnames is None:
            raise SystemExit("inj_metrics.csv has no header")
        alpha_col = _pick_col(r.fieldnames, ["alpha", "alpha_inj"])
        theta_col = _pick_col(r.fieldnames, ["theta", "theta_inj", "theta_data", "theta_value"])
        pp_col = _pick_col(r.fieldnames, ["pp", "PP", "pp_inj", "PP_inj", "pp_value"])
        buf: dict[float, list[tuple[float, float]]] = {}
        for row in r:
            a = float(row[alpha_col])
            buf.setdefault(a, []).append((float(row[theta_col]), float(row[pp_col])))

    out: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for a, pairs in buf.items():
        th = np.asarray([p[0] for p in pairs], dtype=np.float64)
        pp = np.asarray([p[1] for p in pairs], dtype=np.float64)
        out[a] = (th, pp)
    return out


def q_linear(x: np.ndarray, q: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=np.float64), q, method="linear"))


def write_table(
    *,
    path: Path,
    rows: list[dict[str, float]],
) -> None:
    lines: list[str] = []
    lines.append(r"\begin{table}[t]" + "\n")
    lines.append(r"\centering" + "\n")
    lines.append(
        r"\caption{Injection-based sensitivity validation (Section 6). "
        r"The preregistered gate is evaluated at $\alpha=0.10$ and requires CI95 lower bounds $>0$ "
        r"for both $\Delta\theta_{\mathrm{inj}}$ and $\Delta PP_{\mathrm{inj}}$.}" + "\n"
    )
    lines.append(r"\label{tab:inj_effects}" + "\n")
    lines.append(r"\small" + "\n")
    lines.append(r"\setlength{\tabcolsep}{3.5pt}" + "\n")
    lines.append(r"\begin{tabular}{rrrrrrrrr}" + "\n")
    lines.append(r"\toprule" + "\n")
    lines.append(
        r"$\alpha$ & $\theta_{\mathrm{inj}}$ & $\Delta\theta_{\mathrm{inj}}$ & CI95 low & CI95 high & "
        r"$PP_{\mathrm{inj}}$ & $\Delta PP_{\mathrm{inj}}$ & CI95 low & CI95 high \\" + "\n"
    )
    lines.append(r"\midrule" + "\n")
    for x in rows:
        lines.append(
            r"\num{" + f"{x['alpha']:.2f}" + r"} & "
            r"\num{" + f"{x['theta_inj']:.9f}" + r"} & "
            r"\num{" + f"{x['dtheta']:.9f}" + r"} & "
            r"\num{" + f"{x['dtheta_lo']:.9f}" + r"} & "
            r"\num{" + f"{x['dtheta_hi']:.9f}" + r"} & "
            r"\num{" + f"{x['pp_inj']:.6f}" + r"} & "
            r"\num{" + f"{x['dpp']:.6f}" + r"} & "
            r"\num{" + f"{x['dpp_lo']:.6f}" + r"} & "
            r"\num{" + f"{x['dpp_hi']:.6f}" + r"} \\"
            + "\n"
        )
    lines.append(r"\bottomrule" + "\n")
    lines.append(r"\end{tabular}" + "\n")
    lines.append(r"\end{table}" + "\n")
    path.write_text("".join(lines), encoding="utf-8")


def write_section(
    *,
    path: Path,
    run_id: str,
    gate_ok: bool,
) -> None:
    bundle = Path(f"paper/artifacts/inj/{run_id}")
    txt: list[str] = []
    txt.append(
        "We report the preregistered injection-based sensitivity validation (Section 6). "
        "The gate is evaluated at $\\alpha=0.10$ and requires CI95 lower bounds $>0$ for both "
        "$\\Delta\\theta_{\\mathrm{inj}}$ and $\\Delta PP_{\\mathrm{inj}}$ computed against the frozen Gaussian null reference.\n\n"
    )
    if gate_ok:
        txt.append("For this run, the sensitivity gate succeeds at $\\alpha=0.10$.\n\n")
    else:
        txt.append("For this run, the sensitivity gate fails at $\\alpha=0.10$.\n\n")

    txt.append(r"\input{tables/inj_effects_table}" + "\n\n")
    txt.append(r"\paragraph{Artifact traceability.}" + "\n")
    txt.append(
        "All numbers in Table~\\ref{tab:inj_effects} are computed from "
        + r"\path{" + bundle.as_posix() + "/inj_metrics.csv} "
        + "and the frozen null reference stored in "
        + r"\path{" + bundle.as_posix() + "/inj_null_effects.csv}."
        + "\n"
    )
    txt.append(
        "The corresponding run metadata are stored in "
        + r"\path{" + bundle.as_posix() + "/inj_runinfo.json}."
        + "\n"
    )
    txt.append(
        "File digests for the full injection artifact bundle are stored in "
        + r"\path{" + bundle.as_posix() + "/SHA256.txt}."
        + "\n"
    )
    path.write_text("".join(txt), encoding="utf-8")


def ensure_results_includes_inj() -> None:
    p = Path("paper/sections/02_results.tex")
    s = p.read_text(encoding="utf-8")
    if r"\subsection{Injection sensitivity}" in s:
        return
    needle = r"\subsection{HMDIFF negative control}"
    if needle not in s:
        raise SystemExit("paper/sections/02_results.tex does not contain the expected HMDIFF subsection anchor")
    s = s.replace(
        needle,
        r"\subsection{Injection sensitivity}" + "\n" + r"\input{sections/generated/02_results_inj}" + "\n\n" + needle,
    )
    p.write_text(s, encoding="utf-8")
    print("UPDATED: paper/sections/02_results.tex")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    run_id = str(args.run_id).strip()
    bundle = Path(f"paper/artifacts/inj/{run_id}")
    for name in REQ:
        if not (bundle / name).exists():
            raise SystemExit(f"Missing required file: {bundle/name}")

    validate_sha256_ledger(bundle)

    theta_null, pp_null = read_null_arrays(bundle / "inj_null_effects.csv")
    theta_null_ref = float(np.median(theta_null))
    pp_null_ref = float(np.median(pp_null))

    by_alpha = read_inj_by_alpha(bundle / "inj_metrics.csv")

    rows: list[dict[str, float]] = []
    gate_ok = False

    for alpha in sorted(by_alpha.keys()):
        th_arr, pp_arr = by_alpha[alpha]
        theta_inj = float(np.median(th_arr))
        pp_inj = float(np.median(pp_arr))

        dtheta = theta_null_ref - theta_inj
        dtheta_dist = theta_null - theta_inj
        dtheta_lo = q_linear(dtheta_dist, 0.025)
        dtheta_hi = q_linear(dtheta_dist, 0.975)

        dpp = pp_inj - pp_null_ref
        dpp_dist = pp_inj - pp_null
        dpp_lo = q_linear(dpp_dist, 0.025)
        dpp_hi = q_linear(dpp_dist, 0.975)

        if abs(alpha - 0.10) < 1e-12:
            gate_ok = (dtheta_lo > 0.0) and (dpp_lo > 0.0)

        rows.append(
            dict(
                alpha=alpha,
                theta_inj=theta_inj,
                dtheta=dtheta,
                dtheta_lo=dtheta_lo,
                dtheta_hi=dtheta_hi,
                pp_inj=pp_inj,
                dpp=dpp,
                dpp_lo=dpp_lo,
                dpp_hi=dpp_hi,
            )
        )

    Path("paper/tables").mkdir(parents=True, exist_ok=True)
    Path("paper/sections/generated").mkdir(parents=True, exist_ok=True)

    write_table(path=Path("paper/tables/inj_effects_table.tex"), rows=rows)
    print("WROTE: paper/tables/inj_effects_table.tex")

    write_section(
        path=Path("paper/sections/generated/02_results_inj.tex"),
        run_id=run_id,
        gate_ok=gate_ok,
    )
    print("WROTE: paper/sections/generated/02_results_inj.tex")

    ensure_results_includes_inj()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
