# Project Plan v2.9 (Final, Pre-Registration Sealed, Fully Standalone)
## Scale-Dependent Mode Coupling in Planck 2018 Temperature Maps
### Primary Gaussian Null + Fixed End-to-End Attribution Null

---

# 0. Status

This document supersedes v2.8.

Version 2.9 incorporates the final optional robustness clarifications:
- Explicit mask fixation statement.
- Explicit beam smoothing implementation statement.

No statistical, operational, or attributional rule has been modified.
No new degrees of freedom have been introduced.

This document is fully pre-registration sealed and standalone-complete.

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

---

## 2.2 Masking and patch extraction

- Projection: gnomonic.
- Field of view: 12 degrees.
- Patch grid: 256 × 256 pixels.
- Deterministic patch-center list (versioned and frozen).
- Masked pixels are set to 0 before coarse graining.
- Coarse graining uses periodic wrap padding to ensure divisibility.

Mask specification:

The exact mask product is fixed by file name and SHA256 hash and is identical across:
- data,
- primary Gaussian null 𝒩_G,
- injection simulations,
- FFP attribution ensemble 𝒩_FFP.

No alternative masks are permitted.

---

## 2.3 Scale workflow and coarse graining (explicitly frozen)

Fixed scale set:

Score = {1, 2, 4}

For each scale s:

1. Apply periodic wrap padding to the 256 × 256 patch so both dimensions are divisible by s.
2. Apply s × s block-mean averaging.
3. The resulting coarse field has shape (256/s) × (256/s).

This procedure is applied identically to data, 𝒩_G, injections, and 𝒩_FFP.

---

## 2.4 Quantization (fully explicit)

For each patch p and scale s, let x_i denote coarse-grained pixel values.

Define:

xmin = min_i x_i  
xmax = max_i x_i  

Quantization rule:

If xmax = xmin:
- q_i = 0 for all i.

Otherwise:

q_i = floor(255 × (x_i − xmin) / (xmax − xmin))

If floating-point rounding yields 256:
- clamp to 255 deterministically.

No clipping, trimming, percentiles, or global scaling is allowed.

Quantization is strictly patch-internal.

---

## 2.5 Compression and baseline correction

- Byte array constructed in fixed row-major order.
- gzip compression level = 6.
- Backend: Python gzip module (zlib).
- python_version and zlib_version logged.
- Canonical smoke-test array compressed; SHA256 stored.

Define:

bpc(p,s) = compressed_bytes / number_of_cells

Baseline:

bpc0(s) = mean over 3 i.i.d. uniform 8-bit arrays  
(shape identical to q(p,s)).

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

## 2.8 Effect definitions and CI95 (fully explicit)

Primary effect:

Δθ = θ_null_ref − θ_data

Discovery direction:

Δθ > 0

Null reference:

θ_null_ref = median_r θ^(r)

CI95(Δθ) computed from:

Dθ = { θ^(r) − θ_data }_r

CI95(Δθ) = [quantile_0.025(Dθ), quantile_0.975(Dθ)]

---

Auxiliary effect (PP):

Over s ∈ Score:

meanθ = mean_s κ̃θ(s)

PP = 100 × (max_s κ̃θ(s) − min_s κ̃θ(s))  
     / max(|meanθ|, ε)

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

# 3. Primary Null Model 𝒩_G (Decision-Relevant)

Purpose:

Define strictly Gaussian two-point–preserving reference.

---

## 3.1 Cℓ source (hard-fixed)

Planck 2018 best-fit ΛCDM TT theoretical spectrum.

- File path logged.
- SHA256 logged.
- ell_max_used = 6143 fixed.

---

## 3.2 Full-sky synthesis (explicit implementation)

For each realization r:

1. Sample Gaussian a_{ℓm} coefficients consistent with Cℓ^{TT}.
2. Synthesize full-sky map at NSIDE = 2048.
3. Beam smoothing is performed with `healpy.smoothing`
   using a Gaussian beam of FWHM = 10 arcmin in harmonic space.
4. Apply the fixed mask (file name + SHA256 logged).
5. Extract patches via fixed gnomonic projection and center list.
6. Apply identical encoding pipeline.

This ensemble defines all statistical thresholds.

---

# 4. Secondary Attribution Null 𝒩_FFP (Interpretative Only)

Purpose:

Assess whether separation vs 𝒩_G could arise from foreground residuals or instrument structure.

This ensemble does NOT alter discovery thresholds.

---

## 4.1 Simulation source and map type (hard-fixed)

- Planck FFP10 official end-to-end simulation release.
- Component-separated temperature maps labeled SMICA and NILC
  in the official FFP10 distribution.
- Used exactly as provided.
- No custom component separation.
- No reprocessing of frequency maps.
- No alternative SMICA-like variants.

Exact file names and SHA256 hashes are frozen and logged.

---

## 4.2 Ensemble size

N_FFP = 100

The first 100 realizations in canonical release order are used.

No post hoc selection permitted.

---

## 4.3 Processing

Each FFP realization undergoes:

- identical beam smoothing,
- identical mask,
- identical patch extraction,
- identical encoding pipeline,
- identical aggregation and effect computation.

---

## 4.4 Role in inference

𝒩_FFP is interpretative only.

It does not alter decision gates.

---

# 5. Sensitivity Validation (Mandatory and Non-Substitutable)

Injection-based sensitivity validation is binding.

𝒩_FFP does NOT substitute injection validation.

Gaussian confirmation requires sensitivity gate success.

---

# 6. Negative Control

HMDIFF required.

Discovery cannot be declared if HMDIFF satisfies discovery criteria.

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

Every manuscript value must map to artifact path + field name.

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

Version 2.9:

- Fully locks the primary statistical protocol.
- Explicitly defines coarse graining, quantization, and CI distributions.
- Hard-fixes the FFP map type and ensemble size.
- Explicitly fixes mask identity and beam smoothing implementation.
- Maintains strict separation of decision and attribution layers.

No implicit statistical or operational degrees of freedom remain.

This document is pre-registration sealed.