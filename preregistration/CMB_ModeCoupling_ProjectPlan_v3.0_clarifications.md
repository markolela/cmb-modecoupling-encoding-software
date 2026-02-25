# Project Plan v3.0 (Final, Pre-Registration Sealed, Fully Standalone)
## Scale-Dependent Mode Coupling in Planck 2018 Temperature Maps
### Primary Gaussian Null + Fixed End-to-End Attribution Null
### Clarifications Addendum With Frozen External Inputs, Patch Centers, HMDIFF, and Injection Gate

---

# 0. Status

This document supersedes v2.9 by clarifications only.

No statistical, operational, or attributional rule from v2.9 is modified.
No new decision threshold is introduced for discovery.
This document freezes previously underdetermined implementation degrees of freedom that are required for exact reproducibility.

v2.9 remains sealed and unmodified.

This document is pre-registration sealed.

---

# 1. Primary Objective

Test whether Planck 2018 temperature maps (SMICA, NILC) exhibit scale-dependent structure beyond a strictly Gaussian two-point model on the fixed scale set:

Score = {1, 2, 4}

under a fully locked encoding-based diagnostic.

Admissible outcomes:

- Robust separation from the Gaussian two-point null.
- Strong confirmation of Gaussianity, conditional on validated sensitivity.

---

# 2. Binding Primary Statistical System

No modification of any element in Sections 2–7 is allowed.

---

## 2.1 Map and resolution

- Input maps: Planck 2018 SMICA and NILC temperature maps.
- NSIDE = 2048.
- ell_max_used = 6143.

Map field:

- Use I_STOKES only.
- No unit conversion is applied.
- FITS header keywords that describe units are logged in runinfo.

---

## 2.2 Masking and patch extraction

- Projection: gnomonic.
- Field of view: 12 degrees.
- Patch grid: 256 × 256 pixels.
- Deterministic patch-center list is versioned and frozen.
- Masked pixels are set to 0 before coarse graining.
- Coarse graining uses periodic wrap padding to ensure divisibility.

### 2.2.1 Mask specification, product frozen

Mask product is hard-fixed by file name and SHA256, identical across:
- data
- primary Gaussian null 𝒩_G
- injection simulations
- FFP attribution ensemble 𝒩_FFP
- HMDIFF negative control

No alternative masks are permitted.

Frozen mask product:

- File name: COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits
- Local path convention: data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits
- SHA256 must equal the value recorded in:
  preregistration/external_inputs/sha256_mask_planck_pr3_common_int_2048_r3_00.txt

Mask binarization rule:

- A pixel is unmasked if mask_value >= 0.999.
- Otherwise it is masked.
- Masked pixels are set to 0 in the temperature map before patch extraction.

### 2.2.2 Patch center list, versioned and frozen

Frozen patch centers artifact set:

- CSV:
  preregistration/patch_centers/patch_centers_n256_fov12_galactic_v1.csv
- Manifest:
  preregistration/patch_centers/patch_centers_n256_fov12_galactic_v1_manifest.json
- Generator script:
  preregistration/patch_centers/generate_patch_centers_v1.py
- SHA256 ledger:
  preregistration/patch_centers/SHA256.txt

The SHA256 values recorded in preregistration/patch_centers/SHA256.txt are binding.

Coordinate convention:

- Galactic coordinates.
- l_deg in [0, 360).
- b_deg in [-90, 90].
- theta_rad and phi_rad are HEALPix angles in radians.
- theta_rad = colatitude.
- phi_rad = longitude.

Selection algorithm is binding as documented in the manifest JSON and implemented in generate_patch_centers_v1.py.

### 2.2.3 Patch extraction implementation, fully frozen

Input:

- Full sky HEALPix map in RING ordering at NSIDE = 2048.
- Temperature map is already masked to 0 using the frozen mask and threshold rule.

Patch plane grid:

- N = 256.
- FOV_rad = 12 deg converted to radians.
- Pixel spacing:
  delta = FOV_rad / N.
- Pixel center coordinates on tangent plane:
  x_i = (i + 0.5 - N/2) * delta for i in {0,...,N-1}
  y_j = (j + 0.5 - N/2) * delta for j in {0,...,N-1}
- x increases to the east, y increases to the north in the local tangent plane.

Gnomonic inverse projection:

Let the patch center be (theta0, phi0).
For each plane coordinate (x, y), define:
- rho = sqrt(x^2 + y^2)
- c = arctan(rho)

If rho = 0, then theta = theta0 and phi = phi0.

Otherwise:
- sin_theta = cos(c) * cos(theta0) + (y * sin(c) * sin(theta0)) / rho
- theta = arccos(sin_theta_clipped) where sin_theta_clipped is clipped to [-1, 1]
- phi = phi0 + atan2( x * sin(c), rho * sin(theta0) * cos(c) - y * cos(theta0) * sin(c) )

Sampling:

- Use healpy.get_interp_val with nest = False.
- Interpolation is bilinear in (theta, phi).
- Values are returned in float64.

No alternative interpolation orders are permitted.

---

## 2.3 Scale workflow and coarse graining, explicitly frozen

Fixed scale set:

Score = {1, 2, 4}

For each scale s:

1. Apply periodic wrap padding to the 256 × 256 patch so both dimensions are divisible by s.
2. Apply s × s block mean averaging.
3. The resulting coarse field has shape (256/s) × (256/s).

This procedure is applied identically to data, 𝒩_G, injections, 𝒩_FFP, and HMDIFF.

---

## 2.4 Quantization, fully explicit

For each patch p and scale s, let x_i denote coarse-grained pixel values.

Define:

xmin = min_i x_i
xmax = max_i x_i

Quantization rule:

If xmax = xmin:
- q_i = 0 for all i.

Otherwise:
- q_i = floor(255 × (x_i − xmin) / (xmax − xmin))

If floating point rounding yields 256:
- clamp to 255 deterministically.

No clipping, trimming, percentiles, or global scaling is allowed.

Quantization is strictly patch internal.

---

## 2.5 Compression and baseline correction

- Byte array constructed in fixed row major order.
- gzip compression level = 6.
- Backend: Python gzip module and zlib.
- python_version and zlib_version logged.
- Canonical smoke test array compressed and SHA256 stored.

Define:

bpc(p,s) = compressed_bytes / number_of_cells

Baseline:

bpc0(s) = mean over 3 i.i.d. uniform 8-bit arrays
shape identical to q(p,s)

RNG:

- numpy.random.Generator(PCG64)
- baseline seed fixed and logged.

Encoding proxy:

κ(p,s) = (bpc(p,s) − bpc0(s)) × ln(256)

All logarithms are natural logarithms.

---

## 2.6 Aggregation

Median across patches:

κ̃(s) = median_p κ(p,s)

---

## 2.7 Trend model and detrending

Trend model:

m(s; a, θ) = a + θ log(s)

Fit:

- Unweighted OLS
- Fit domain: s ∈ Score
- No alternative model families
- No robust regression
- No weighting

Let θ_data be the fitted slope.

Detrended curve:

κ̃θ(s) = κ̃(s) − θ_data · log(s)

---

## 2.8 Effect definitions and CI95, fully explicit

Primary effect:

Δθ = θ_null_ref − θ_data

Discovery direction:

Δθ > 0

Null reference:

θ_null_ref = median_r θ^(r)

CI95(Δθ) computed from:

Dθ = { θ^(r) − θ_data }_r

CI95(Δθ) = [quantile_0.025(Dθ), quantile_0.975(Dθ)]

Auxiliary effect (PP):

Over s ∈ Score:

meanθ = mean_s κ̃θ(s)

PP = 100 × (max_s κ̃θ(s) − min_s κ̃θ(s)) / max(|meanθ|, ε)

ε = 1e−12

For each null draw r:

PP^(r) computed identically.

Null reference:

PP_null_ref = median_r PP^(r)

Define:

ΔPP = PP_data − PP_null_ref

CI95(ΔPP) from:

DPP = { PP_data − PP^(r) }_r

CI95(ΔPP) = [quantile_0.025(DPP), quantile_0.975(DPP)]

---

## 2.9 Null budget

Primary Gaussian null:

N_null = 512

---

# 3. Primary Null Model 𝒩_G (Decision Relevant)

Purpose:

Define strictly Gaussian two-point preserving reference.

---

## 3.1 Cℓ source, hard fixed

Planck 2018 best fit ΛCDM TT theoretical spectrum.

Frozen file name:

- COM_PowerSpect_CMB-base-plikHM_TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt

Local path convention:

- data/external/planck_pr3/COM_PowerSpect_CMB-base-plikHM_TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt

SHA256 binding rule:

- SHA256 must be recorded in preregistration/external_inputs/ as a dedicated sha256 text file.
- The value recorded there is binding.

ell_max_used = 6143 is fixed.

---

## 3.2 Full sky synthesis, explicit implementation frozen

For each realization r:

1. Sample Gaussian a_lm coefficients consistent with C_l^{TT} up to ell_max_used.
2. Synthesize full sky map at NSIDE = 2048 in RING ordering.
3. Beam smoothing is performed with healpy.smoothing using a Gaussian beam of FWHM = 10 arcmin in harmonic space.
4. Apply the fixed mask using the threshold rule in Section 2.2.1.
5. Extract patches via the frozen gnomonic projection in Section 2.2.3 and the frozen center list in Section 2.2.2.
6. Apply the identical encoding pipeline.

Smoothing parameters are frozen:

- fwhm = 10 arcmin converted to radians.
- lmax = ell_max_used.
- pol = False.
- nest = False.

All map arrays are processed in float64.

RNG and determinism for null draws are frozen:

- A per realization integer seed is derived as seed_r = null_seed + r, with r starting at 0.
- The seed_r value is logged for every realization.
- The RNG mechanism used by the synthesis implementation must be fully logged in runinfo.

This ensemble defines all statistical thresholds.

---

# 4. Secondary Attribution Null 𝒩_FFP (Interpretative Only)

Purpose:

Assess whether separation vs 𝒩_G could arise from foreground residuals or instrument structure.

This ensemble does not alter discovery thresholds.

---

## 4.1 Simulation source and map type, hard fixed

- Planck FFP10 official end to end simulation release.
- Component separated temperature maps labeled SMICA and NILC in the official FFP10 distribution.
- Used exactly as provided.
- No custom component separation.
- No reprocessing of frequency maps.
- No alternative SMICA like variants.

FFP file freezing rule:

- Exact file names and SHA256 hashes are frozen in a committed manifest file:
  preregistration/external_inputs/ffp10_manifest_v1.json
- This manifest must list:
  realization index, file name, local path, SHA256, and map field used.

---

## 4.2 Ensemble size

N_FFP = 100

The first 100 realizations in canonical release order are used.

No post hoc selection permitted.

---

## 4.3 Processing

Each FFP realization undergoes:

- identical beam smoothing
- identical mask
- identical patch extraction
- identical encoding pipeline
- identical aggregation and effect computation

---

## 4.4 Role in inference

𝒩_FFP is interpretative only.

It does not alter decision gates.

---

# 5. Sensitivity Validation (Mandatory and Non Substitutable)

Injection based sensitivity validation is binding.

𝒩_FFP does not substitute injection validation.

Gaussian confirmation requires sensitivity gate success.

This section freezes the injection suite that was referenced in v2.9 but not operationally specified.

---

## 5.1 Injection suite, frozen

All injections use the same pipeline as 𝒩_G and differ only by a deterministic modification applied before masking.

Injection type A, large scale modulation that induces mode coupling:

1. Draw a base Gaussian full sky map T_G via the 𝒩_G synthesis.
2. Draw an independent low ell modulation field M:
   - M is synthesized from the same C_l^{TT} file after zeroing all multipoles with ell > 20.
   - M is smoothed with the same 10 arcmin beam.
   - M is normalized to mean 0 and standard deviation 1 over unmasked pixels.
3. Form the injected map:
   T_inj = (1 + alpha * M) * T_G

Alpha grid is frozen:

- alpha in {0.05, 0.10, 0.15}

Ensemble size per alpha is frozen:

- N_inj = 128

---

## 5.2 Sensitivity gate rule, frozen

The gate is evaluated at alpha = 0.10.

Compute Δθ_inj and ΔPP_inj versus the same 𝒩_G reference definition as in Sections 2.8 and 2.9.

Gate success requires both conditions:

1. CI95(Δθ_inj) lower bound is greater than 0.
2. CI95(ΔPP_inj) lower bound is greater than 0.

If the gate fails, the study may still report separation claims, but it may not claim strong confirmation of Gaussianity in case of null like outcomes.

All injection artifacts and gate outputs are written and hashed under the artifact policy.

---

# 6. Negative Control

HMDIFF is required.

Discovery cannot be declared if HMDIFF satisfies discovery criteria.

This section freezes the HMDIFF construction that was referenced in v2.9 but not operationally specified.

---

## 6.1 HMDIFF construction, frozen

For each map product family separately:

- SMICA:
  - Use COM_CMB_IQU-smica_2048_R3.00_hm1.fits and COM_CMB_IQU-smica_2048_R3.00_hm2.fits
- NILC:
  - Use COM_CMB_IQU-nilc_2048_R3.00_hm1.fits and COM_CMB_IQU-nilc_2048_R3.00_hm2.fits

Use I_STOKES field only.

Define the half mission difference map:

- HMDIFF = 0.5 * (HM1 - HM2)

Then apply:

- identical beam smoothing
- identical mask
- identical patch extraction
- identical encoding pipeline
- identical aggregation and effect computation

File names and SHA256 hashes are frozen and logged under preregistration/external_inputs/ using the same binding rule as for the mask.

---

# 7. Artifact Policy

Primary ensemble artifacts:

- metrics.csv
- metrics_per_patch.csv
- null_effects.csv
- summary.json
- null_summary.json
- runinfo.json
- gzip_smoketest_hash.txt

FFP artifacts:

- ffp_metrics.csv
- ffp_summary.json
- ffp_runinfo.json

Injection artifacts:

- inj_summary.json
- inj_null_summary.json
- inj_runinfo.json
- inj_metrics.csv
- inj_null_effects.csv

HMDIFF artifacts:

- hmdiff_summary.json
- hmdiff_runinfo.json
- hmdiff_metrics.csv

Every manuscript value must map to artifact path and field name.

Runinfo must include:

- git_commit
- python_version
- numpy_version
- scipy_version
- healpy_version
- astropy_version
- zlib_version
- platform string
- all config values
- all seeds and derived per realization seeds
- all external input file paths and SHA256 hashes
- patch center list SHA256 and manifest SHA256

---

# 8. Scope Discipline

Restricted to temperature maps.

Excluded:

- polarization
- cosmological parameter inference
- bispectrum fitting
- ontological interpretation without statistical separation

---

# 9. Final Statement

Version 3.0:

- Fully locks all v2.9 decision relevant statistics and definitions.
- Freezes the concrete mask product and thresholding rule.
- Freezes the patch center list artifacts and the exact projection and interpolation.
- Freezes HMDIFF construction from half mission splits.
- Freezes the injection sensitivity suite and the gate criterion.
- Strengthens auditability by binding external input hashes via committed ledgers.

This document is pre-registration sealed.