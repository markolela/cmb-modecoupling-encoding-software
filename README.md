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
