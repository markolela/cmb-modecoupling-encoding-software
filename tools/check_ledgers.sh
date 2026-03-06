#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "REPO_ROOT: ${repo_root}"
echo

echo "== preregistration patch_centers SHA256 =="
sha256sum -c preregistration/patch_centers/SHA256.txt
echo

echo "== preregistration SHA256, relative paths =="
(
  cd preregistration
  sha256sum -c SHA256.txt
)
echo

echo "== external_inputs env lock ledger =="
got="$(sha256sum preregistration/external_inputs/requirements_lock_wsl_ubuntu.txt | awk '{print $1}')"
exp="$(tr -d '\r\n ' < preregistration/external_inputs/sha256_env_lock_wsl_ubuntu.txt)"
echo "GOT: ${got}"
echo "EXP: ${exp}"
test "${got}" = "${exp}" && echo "OK" || { echo "MISMATCH"; exit 2; }
echo

echo "== external_inputs FFP manifest ledger =="
if [ -f preregistration/external_inputs/ffp10_manifest_v1.json ] && [ -f preregistration/external_inputs/sha256_ffp10_manifest_v1.txt ]; then
  got="$(sha256sum preregistration/external_inputs/ffp10_manifest_v1.json | awk '{print $1}')"
  exp="$(tr -d '\r\n ' < preregistration/external_inputs/sha256_ffp10_manifest_v1.txt)"
  echo "GOT: ${got}"
  echo "EXP: ${exp}"
  test "${got}" = "${exp}" && echo "OK" || { echo "MISMATCH"; exit 2; }
else
  echo "SKIP, manifest oder ledger fehlt"
fi
echo

echo "ALL LEDGER CHECKS PASSED"
