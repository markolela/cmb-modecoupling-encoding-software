# Scripts

## Planck beam harmonization (10 arcmin)

The paper uses beam harmonization to 10 arcmin FWHM prior to patch extraction.
The following scripts are used for this preprocessing step.

### Generic harmonization (beam + NSIDE)

Example:
python scripts/harmonize_beam_nside.py \
  --in <input_fits_1> <input_fits_2> \
  --fwhm-arcmin 10 \
  --nside <nside> \
  --outdir data/raw/astro/planck/harmonized

### Half-mission harmonization helper (SMICA HM1/HM2)

Example:
python scripts/harmonize_planck_hm_sm10am.py \
  --in-hm1 <hm1_fits> \
  --in-hm2 <hm2_fits> \
  --outdir data/raw/astro/planck/harmonized
