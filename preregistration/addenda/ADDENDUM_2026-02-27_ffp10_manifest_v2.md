# Addendum 2026-02-27: FFP10 manifest v2 snapshot

## What was missing or defective
The repository contained only a placeholder FFP10 manifest v1 with an empty entry list.
The local FFP10 data root initially contained no FFP10 FITS inputs.

## Why execution was blocked
Section 5 requires a manifest snapshot listing N_FFP = 100 local inputs with SHA256.
With no entries, Section 5 could not be executed.

## Exact fix applied and why it restores the frozen protocol
We add a new manifest version:
- preregistration/external_inputs/ffp10_manifest_v2.json
and its SHA256 ledger:
- preregistration/external_inputs/sha256_ffp10_manifest_v2.txt

This does not modify any decision relevant section. It only supplies missing secondary interpretative inputs required for Section 5.

Selection rule for v2.
We enumerate exactly 100 FITS files located directly under data/external/ffp10/ in lexicographic filename order.
Each entry records local_path, SHA256, and field=0.

## Exact inputs used
All inputs are under data/external/ffp10/ and are enumerated with SHA256 in ffp10_manifest_v2.json.

## Output directory for the backfilled execution
Section 5 outputs must be written to a new dedicated Section 5 output directory and must cite the manifest snapshot pair:
ffp_manifest_used.json and ffp_manifest_ledger_sha256.txt.
