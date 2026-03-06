#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "REPO_ROOT: ${repo_root}"
echo

echo "== git status =="
git status
echo

echo "== ledger checks =="
./tools/check_ledgers.sh
echo

echo "== presence of primary external inputs =="
for f in \
  data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_full.fits \
  data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_full.fits \
  data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits \
  data/external/planck_pr3/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt \
  data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_hm1.fits \
  data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_hm2.fits \
  data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_hm1.fits \
  data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_hm2.fits \
  gzip_smoketest_hash.txt
do
  if [ -f "$f" ]; then
    echo "OK: $f"
  else
    echo "MISSING: $f"
    exit 2
  fi
done

echo
echo "HEALTH CHECK PASSED"
