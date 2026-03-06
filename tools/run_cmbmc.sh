#!/usr/bin/env bash
set -u

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}" || { echo "Repo root not found"; pwd; exit 0; }

# Use venv python, and make src importable.
PY="${repo_root}/.venv/bin/python"
if [ ! -x "${PY}" ]; then
  echo "Missing venv python at ${PY}"
  exit 0
fi

PYTHONPATH="${repo_root}/src" "${PY}" -m cmbmc.cli "$@"
