# cmb-modecoupling-encoding-software

Preregistered, artifact-traceable pipeline for scale-dependent mode coupling diagnostics in Planck 2018 temperature maps (SMICA/NILC).

This repository contains:
- src/ (cmbmc package)
- tools/ (runners and paper update scripts)
- preregistration/ (sealed protocol, addenda, and external input ledgers)

Paper sources and large artifacts are released separately (arXiv / Zenodo).

Quick start:
1) Create venv: python -m venv .venv && source .venv/bin/activate
2) Install: pip install -r preregistration/external_inputs/requirements_lock_wsl_ubuntu.txt
3) Run: ./tools/run_cmbmc.sh primary --out /tmp/cmbmc_primary_out

## Citation and DOIs

### Software (this repository)
Zenodo DOI (v0.1.0): 10.5281/zenodo.18883337

Use the GitHub citation panel or CITATION.cff for metadata.

### Paper B artifacts (run bundles)
Zenodo DOI (version): 10.5281/zenodo.18883260
Zenodo DOI (all versions): 10.5281/zenodo.18883259

This dataset contains the archived Paper B artifact bundles for the registered runs (primary, ffp, hmdiff, inj), including SHA256 ledgers and paper/runs.json.
