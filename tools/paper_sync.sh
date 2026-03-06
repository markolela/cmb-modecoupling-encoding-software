#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

PY=".venv/bin/python"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

runs_json="paper/runs.json"
test -f "$runs_json" || { echo "FEHLT: $runs_json"; exit 1; }

primary="$($PY - <<'PY'
import json
print(json.load(open("paper/runs.json","r",encoding="utf-8"))["primary"])
PY
)"
ffp="$($PY - <<'PY'
import json
print(json.load(open("paper/runs.json","r",encoding="utf-8"))["ffp"])
PY
)"
hmdiff="$($PY - <<'PY'
import json
print(json.load(open("paper/runs.json","r",encoding="utf-8"))["hmdiff"])
PY
)"
inj="$($PY - <<'PY'
import json
print(json.load(open("paper/runs.json","r",encoding="utf-8"))["inj"])
PY
)"

echo "RUNS:"
echo "  primary: $primary"
echo "  ffp    : $ffp"
echo "  hmdiff : $hmdiff"
echo "  inj    : $inj"
echo

# Vorchecks: Artefakte müssen existieren
test -f "paper/artifacts/primary/${primary}/summary.json" || { echo "FEHLT: paper/artifacts/primary/${primary}/summary.json"; exit 1; }
test -f "paper/artifacts/ffp/${ffp}/ffp_summary.json" || { echo "FEHLT: paper/artifacts/ffp/${ffp}/ffp_summary.json"; exit 1; }
test -f "paper/artifacts/hmdiff/${hmdiff}/hmdiff_summary.json" || { echo "FEHLT: paper/artifacts/hmdiff/${hmdiff}/hmdiff_summary.json"; exit 1; }
test -f "paper/artifacts/inj/${inj}/inj_metrics.csv" || { echo "FEHLT: paper/artifacts/inj/${inj}/inj_metrics.csv"; exit 1; }

# Generated Sections + Tabellen neu erzeugen
./tools/paper_update_primary_results.py --run-id "$primary"
./tools/paper_update_ffp_results.py --run-id "$ffp" --primary-run-id "$primary"
./tools/paper_update_hmdiff_results.py --run-id "$hmdiff"
./tools/paper_update_inj_results.py --run-id "$inj"


# 02_results.tex kanonisch setzen (verhindert doppelte Includes durch einzelne Update-Skripte)
cat > paper/sections/02_results.tex <<'TEX'
\section{Results}

\subsection{Primary results}
\input{sections/generated/02_results_primary}

\subsection{FFP10 attribution ensemble}
\input{sections/generated/02_results_ffp}

\subsection{Injection-based sensitivity validation}
\input{sections/generated/02_results_inj}

\subsection{HMDIFF negative control}
\input{sections/generated/02_results_hmdiff}
TEX

# LaTeX build artifacts löschen (verhindert stale duplicate-label warnings)
rm -f paper/main.aux paper/main.out paper/main.log paper/main.toc paper/main.pdf

# PDF 2x bauen (Referenzen)
( cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex )
( cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex )

# Harte Fehler prüfen
grep -nE "Fatal error|Missing \\$ inserted|Undefined control sequence|LaTeX Error" paper/main.log && {
  echo "HARTE LATEX FEHLER GEFUNDEN. Siehe paper/main.log."
  exit 2
} || true

test -f paper/main.pdf || { echo "FEHLT: paper/main.pdf"; exit 1; }

echo "OK: paper sync + build erfolgreich"
ls -la paper/main.pdf
