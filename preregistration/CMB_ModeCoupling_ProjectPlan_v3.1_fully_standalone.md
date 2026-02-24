# File: preregistration/CMB_ModeCoupling_ProjectPlan_v3.1_fully_standalone.md
# Project Plan v3.1 (Final, Pre-Registration Sealed, Fully Standalone)
## Scale-Dependent Mode Coupling in Planck 2018 Temperature Maps
### Primary Gaussian Null + Fixed End-to-End Attribution Null
### Fully Frozen Inputs, Seeds, RNG, Index Conventions, Decision Rules, and Injection Suite

---

# 0. Status

This document supersedes v3.0 by clarification and freezing only.

v2.9 and v3.0 remain sealed and unmodified.

No new reported quantity is introduced.
No effect definition is changed.
No thresholds are loosened.
This document freezes previously underdetermined degrees of freedom that are required for exact reproduction from the plan alone.

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

No modification of any element in Sections 2–8 is allowed.

---

## 2.0 Frozen global constants and seeds

All seeds are binding integers.

- baseline_seed = 2026022401
- null_seed = 2026022402
- inj_seed = 2026022403

All runs must log the seeds and derived per-realization seeds in runinfo.

---

## 2.1 Data products (hard-fixed)

The exact data products are frozen by file name, local path convention, and SHA256.

SMICA full mission:

- File name: COM_CMB_IQU-smica_2048_R3.00_full.fits
- Local path convention:
  data/external/planck_pr3/COM_CMB_IQU-smica_2048_R3.00_full.fits
- SHA256 must equal the value recorded in:
  preregistration/external_inputs/sha256_data_smica_pr3_full.txt

NILC full mission:

- File name: COM_CMB_IQU-nilc_2048_R3.00_full.fits
- Local path convention:
  data/external/planck_pr3/COM_CMB_IQU-nilc_2048_R3.00_full.fits
- SHA256 must equal the value recorded in:
  preregistration/external_inputs/sha256_data_nilc_pr3_full.txt

Map field (binding):

- Use I_STOKES only.
- I_STOKES is loaded as FITS field index 0.
- The binding load call is:
  hp.read_map(path, field=0, nest=False, dtype=numpy.float64, verbose=False)

No unit conversion is applied.

Resolution:

- NSIDE = 2048.
- ell_max_used = 6143.

---

## 2.2 Masking and beam smoothing (data pipeline is binding)

### 2.2.1 Mask specification (hard-fixed)

Mask product is hard-fixed by file name and SHA256 and is identical across:
- data,
- primary Gaussian null 𝒩_G,
- injection simulations,
- FFP attribution ensemble 𝒩_FFP,
- HMDIFF negative controls.

Frozen mask product:

- File name: COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits
- Local path convention:
  data/external/planck_pr3/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits
- SHA256 must equal the value recorded in:
  preregistration/external_inputs/sha256_mask_planck_pr3_common_int_2048_r3_00.txt

Binding load call (mask):

- The mask map is loaded as:
  hp.read_map(path, field=0, nest=False, dtype=numpy.float64, verbose=False)

Mask binarization rule (binding):

- A pixel is unmasked if mask_value >= 0.999.
- Otherwise it is masked.

Mask application rule (binding):

- Masked pixels are set to 0 in the temperature map.
- This is done immediately before patch extraction.

### 2.2.2 Beam smoothing (binding for all pipelines)

Beam smoothing is mandatory for:
- data maps,
- 𝒩_G realizations,
- injections,
- FFP attribution realizations,
- HMDIFF maps.

Smoothing is performed with healpy.smoothing using a Gaussian beam of FWHM = 10 arcmin.

Frozen smoothing parameters and binding call:

- Define fwhm_10arcmin_rad = 10 arcmin converted to radians.

The binding smoothing call is:

hp.smoothing(
  map_in,
  fwhm=fwhm_10arcmin_rad,
  lmax=ell_max_used,
  pol=False,
  nest=False,
  verbose=False
)

Order of operations is binding:

1. Load map (data, null, injection, FFP, or HM split).
2. Apply healpy.smoothing with the frozen parameters above.
3. Apply the frozen mask to 0 using the rule in 2.2.1.
4. Extract patches via the frozen patch extraction in 2.3.

No alternative ordering is permitted.

---

## 2.3 Patch centers (versioned and frozen)

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
- theta_rad is colatitude.
- phi_rad is longitude.

---

## 2.4 Patch extraction implementation (fully frozen, index-level)

Projection:

- Gnomonic.
- Field of view: 12 degrees.
- Patch grid: 256 × 256 pixels.

Array convention is binding:

- Patch array is A[j, i].
- j indexes y (north direction), i indexes x (east direction).
- A.shape = (256, 256).
- Row-major byte stream is constructed from A.ravel(order="C"),
  which iterates i fastest within each row j.

Plane grid definition (binding):

Let N = 256 and FOV_rad = 12 deg converted to radians.
Pixel spacing:

- delta = FOV_rad / N

Pixel center coordinates:

- x_i = (i + 0.5 - N/2) * delta for i in {0,...,N-1}
- y_j = (j + 0.5 - N/2) * delta for j in {0,...,N-1}

Gnomonic inverse projection is binding:

Let the patch center be (theta0, phi0).
For each (x, y):

- rho = sqrt(x^2 + y^2)
- c = arctan(rho)

If rho = 0:
- theta = theta0
- phi = phi0

Else:
- cos_theta = cos(c) * cos(theta0) + (y * sin(c) * sin(theta0)) / rho
- theta = arccos(clip(cos_theta, -1, 1))
- phi = phi0 + atan2( x * sin(c), rho * sin(theta0) * cos(c) - y * cos(theta0) * sin(c) )

Phi wrap rule (binding):

- Before interpolation, enforce:
  phi = phi mod (2*pi)

Sampling rule (binding):

- Use healpy.get_interp_val with nest = False.
- Bilinear interpolation in (theta, phi).
- Output dtype is float64.

---

## 2.5 Scale workflow and coarse graining (explicitly frozen, index-level)

Fixed scale set:

Score = {1, 2, 4}

For each scale s:

1. If 256 mod s != 0, apply periodic wrap padding to make both dimensions divisible by s.
   For the frozen Score set this is a no-op, but the rule is binding.
2. Apply s × s block-mean averaging with fixed block segmentation:
   blocks start at A[0,0] and proceed contiguously in both axes.
3. Axis convention is binding:
   axis 0 is y (rows), axis 1 is x (columns).
4. Block mean is defined by:
   A_s = A.reshape(256/s, s, 256/s, s).mean(axis=(1,3))
5. The resulting coarse field has shape (256/s, 256/s).

This procedure is applied identically to data, 𝒩_G, injections, 𝒩_FFP, and HMDIFF.

---

## 2.6 Quantization (fully explicit)

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

Byte order and flattening is binding:

- Construct byte array from q (uint8) using q.ravel(order="C"),
  where q is indexed as q[j,i] with the same A[j,i] convention.

---

## 2.7 Compression and baseline correction (fully frozen)

Compression is fully deterministic.

Gzip rule is binding:

- Compress exclusively with:
  gzip.compress(data_bytes, compresslevel=6, mtime=0)
- Define:
  gzip_bytes = gzip.compress(data_bytes, compresslevel=6, mtime=0)
- Then:
  compressed_bytes = len(gzip_bytes)
- The value compressed_bytes is:
  len(gzip_bytes)
  including gzip header and footer.

Byte stream rule is binding:

- For each patch p and scale s, data_bytes is the raw bytes of q as defined in Section 2.6,
  constructed from q.ravel(order="C") where q is indexed q[j,i].

Define number_of_cells (binding):

- For scale s:
  H = 256/s
  W = 256/s
  number_of_cells = H*W

Then:

bpc(p,s) = compressed_bytes / number_of_cells

Baseline bpc0(s) is frozen and computed once per scale, not per patch.

Scale order is binding:

- Iterate scales in this exact order:
  s_list = [1, 2, 4]

Baseline construction is binding:

- RNG: numpy.random.Generator(numpy.random.PCG64(baseline_seed))
- For each scale s in s_list:
  - Let H = 256/s and W = 256/s.
  - Set number_of_cells = H*W.
  - Draw exactly three arrays in a fixed loop:
    U_k = rng.integers(0, 256, size=(H, W), dtype=numpy.uint8) for k in {1,2,3}
  - For each U_k:
    - data_bytes = U_k.ravel(order="C").tobytes()
    - gzip_bytes = gzip.compress(data_bytes, compresslevel=6, mtime=0)
    - bpc_k(s) = len(gzip_bytes) / number_of_cells
  - Set:
    bpc0(s) = mean_k bpc_k(s)

Encoding proxy:

κ(p,s) = (bpc(p,s) − bpc0(s)) × ln(256)

All logarithms are natural logarithms.

Canonical gzip smoke test is binding.

Smoke test byte sequence:

- Define an array of 1024 bytes:
  smoke = [0,1,2,...,255] repeated exactly 4 times.
- Define:
  smoke_bytes = numpy.asarray(smoke, dtype=numpy.uint8).tobytes()

Smoke test procedure:

1. gzip_bytes = gzip.compress(smoke_bytes, compresslevel=6, mtime=0)
2. smoke_sha256 = SHA256(gzip_bytes)

The value smoke_sha256 is stored as a text artifact:

- gzip_smoketest_hash.txt

---

## 2.8 Aggregation (fully frozen)

Scale order is binding:

- s_list = [1, 2, 4]

Let K be the patch by scale matrix:

- K[p, t] = κ(p, s_list[t]) as float64

Median aggregation is binding:

- κ̃ = numpy.median(K, axis=0)

This yields κ̃(t) for t in {0,1,2} corresponding to s_list.

---

## 2.9 Trend model and detrending (OLS call frozen)

Scale order is binding:

- s_list = [1, 2, 4]

Define:

- y = κ̃ as float64 vector of length 3
- x = numpy.log(numpy.asarray(s_list, dtype=numpy.float64))

Design matrix is binding:

- X is a float64 matrix of shape (3, 2)
- X[:,0] = 1
- X[:,1] = x

OLS routine is binding:

- beta = numpy.linalg.lstsq(X, y, rcond=None)[0]
- a = beta[0]
- θ_data = beta[1]

Detrended curve is binding:

- κ̃θ(s_list[t]) = κ̃[t] − θ_data * log(s_list[t])

---

## 2.10 Effect definitions and CI95 (fully explicit, quantiles frozen)

Primary effect:

Δθ = θ_null_ref − θ_data

Discovery direction:

Δθ > 0

Null reference:

θ_null_ref = median_r θ^(r)

CI95(Δθ) computed from:

Dθ = { θ^(r) − θ_data }_r

Quantile operator is binding:

- quantile_q(x) = numpy.quantile(x, q, method="linear")

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

## 2.11 Null budget

Primary Gaussian null:

N_null = 512

---

# 3. Decision rules (explicitly frozen)

The study has two primary data products: SMICA and NILC.

The decision claim "Robust separation from the Gaussian two-point null" is declared if and only if:

1. For SMICA: CI95(Δθ)_lower > 0
2. For NILC:  CI95(Δθ)_lower > 0

This is a conjunctive rule. Both products must satisfy the criterion.

ΔPP is auxiliary and is reported with CI95, but it is not a decision gate.

Strong confirmation of Gaussianity is permitted only if the sensitivity gate in Section 6 succeeds.

---

# 4. Primary Null Model 𝒩_G (Decision-Relevant)

Purpose:

Define strictly Gaussian two-point–preserving reference.

---

## 4.1 Cℓ source (hard-fixed)

Planck 2018 best-fit ΛCDM TT theoretical spectrum.

Frozen file name:

- COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt

Local path convention:

- data/external/planck_pr3/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt

SHA256 must equal the value recorded in:

- preregistration/external_inputs/sha256_cl_theory_planck_pr3_minimum_theory_r3_01.txt

ell_max_used = 6143 fixed.

Cℓ parsing rule (binding):

- Read the file as ASCII text, line by line.
- Ignore empty lines.
- Ignore any line where, after strip, the first character is not a digit.
- Ignore any line where, after strip, the first character is '#'.

For each remaining line:

- Parse whitespace-separated tokens as floats.
- Interpret the first numeric token as ell.
- Interpret the second numeric token as D_ell^TT in microK^2.

Construct the array cl_tt as follows:

- Initialize cl_tt with length ell_max_used + 1 as float64 zeros.
- For each parsed (ell, D_ell^TT):
  - If ell < 2, skip.
  - If ell > ell_max_used, skip.
  - Convert to C_ell^TT via:
    C_ell^TT = D_ell^TT * 2*pi / (ell*(ell+1))
  - Set cl_tt[ell] = C_ell^TT
- Set cl_tt[0] = 0 and cl_tt[1] = 0.

This yields the cl_tt array used by the synthesis.

---

## 4.2 Full-sky synthesis (RNG and function call frozen)

Synthesis implementation is binding.

A single Gaussian null ensemble is generated once and is reused for both data products (SMICA and NILC).
The identical set {θ^(r)} defines θ_null_ref and Dθ for both products.
No separate null ensembles per product are permitted.

For each realization r in {0,...,N_null-1}:

1. Define seed_r = null_seed + r
2. Set the legacy NumPy global RNG seed:
   numpy.random.seed(seed_r)
3. Determinism rule (binding):
   - After numpy.random.seed(seed_r), and until hp.synfast returns,
     no other RNG calls are permitted in the process.

4. Synthesize a full-sky Gaussian map using healpy.synfast with fully frozen parameters:

   hp.synfast(
     cl_tt,
     nside=2048,
     lmax=ell_max_used,
     pol=False,
     new=True,
     verbose=False,
     pixwin=False,
     fwhm=0.0
   )

5. Apply beam smoothing with hp.smoothing using the binding call in Section 2.2.2.
6. Apply the fixed mask-to-0 step using Section 2.2.1.
7. Extract patches via Sections 2.3 and 2.4.
8. Apply the identical encoding pipeline.

All arrays are processed in float64.

This ensemble defines all statistical thresholds.

---

# 5. Secondary Attribution Null 𝒩_FFP (Interpretative Only)

Purpose:

Assess whether separation vs 𝒩_G could arise from foreground residuals or instrument structure.

This ensemble does not alter decision gates.

FFP manifest freezing rule (binding):

- Exact file names and SHA256 hashes are frozen in:
  preregistration/external_inputs/ffp10_manifest_v1.json

The manifest must list:

- realization index
- file name
- local path
- SHA256
- map field used

Ensemble size:

- N_FFP = 100
- The first 100 realizations in canonical release order are used.

Processing:

Each FFP realization undergoes the same order:

Binding load call (FFP realizations):

- Each FFP map is loaded as:
  hp.read_map(path, field=0, nest=False, dtype=numpy.float64, verbose=False)

1. Load map
2. Smooth
3. Mask to 0
4. Patch
5. Encode
6. Aggregate and compute effects

---

# 6. Sensitivity Validation (Injection Suite and Gate Frozen)

Injection-based sensitivity validation is binding.

𝒩_FFP does not substitute injection validation.

Gaussian confirmation requires gate success.

This section freezes injection definitions, normalization, and seeding.

---

## 6.1 Injection suite (frozen)

All injections use the same pipeline as 𝒩_G and differ only by a deterministic modification applied before the mask-to-0 step.

Injection type A, multiplicative low-ell modulation:

Alpha indexing (binding):

- alpha_list = [0.05, 0.10, 0.15]
- alpha_index is the index of alpha in alpha_list, in this exact order.
- Therefore:
  alpha_index(0.05) = 0
  alpha_index(0.10) = 1
  alpha_index(0.15) = 2

All uses of alpha_index in seeding are with this definition.

1. Draw a base Gaussian full-sky map T_G as synfast plus smoothing only.

   Binding seeding:
   - seed_T = inj_seed + 100000*alpha_index + r
   - numpy.random.seed(seed_T)

   Determinism rule (binding):
   - After numpy.random.seed(seed_T), and until hp.synfast returns,
     no other RNG calls are permitted in the process.

   Binding synthesis:
   - T_G = hp.synfast(
       cl_tt,
       nside=2048,
       lmax=ell_max_used,
       pol=False,
       new=True,
       verbose=False,
       pixwin=False,
       fwhm=0.0
     )

   No smoothing is applied to T_G at this stage.
   The only mandatory beam smoothing step for injections is applied to T_inj after it is formed.
   No mask-to-0 is applied to T_G at this stage.

2. Draw an independent low-ell modulation field M as synfast plus smoothing only.

   Binding low-ell spectrum:
   - Construct cl_tt_low = cl_tt, but set cl_tt_low[ell] = 0 for ell > 20.

   Binding seeding:
   - seed_M = inj_seed + 500000 + 100000*alpha_index + r
   - numpy.random.seed(seed_M)

   Determinism rule (binding):
   - After numpy.random.seed(seed_M), and until hp.synfast returns,
     no other RNG calls are permitted in the process.

   Binding synthesis:
   - M = hp.synfast(
       cl_tt_low,
       nside=2048,
       lmax=ell_max_used,
       pol=False,
       new=True,
       verbose=False,
       pixwin=False,
       fwhm=0.0
     )

   No smoothing is applied to M at this stage.
   The only mandatory beam smoothing step for injections is applied to T_inj after it is formed.
   No mask-to-0 is applied to M at this stage.

3. Define the unmasked pixel set U using the frozen mask threshold (mask_value >= 0.999).
   The mask is used only as a selector for U.
   No pixels are set to 0 at this step.

4. Normalize M over U:

   - mu = mean_{pix in U} M[pix]
   - sigma = std_{pix in U} M[pix]
   - M_norm = (M - mu) / sigma

5. Form the injected map:

   T_inj = (1 + alpha * M_norm) * T_G

6. Apply beam smoothing to T_inj using the binding call in Section 2.2.2.

7. Apply mask-to-0 to T_inj using Section 2.2.1.
8. Extract patches and run the identical encoding pipeline.

Alpha grid (binding):

- alpha in {0.05, 0.10, 0.15}

Ensemble size per alpha (binding):

- N_inj = 128

---

## 6.2 Sensitivity gate rule (binding)

The gate is evaluated at alpha = 0.10.

Compute Δθ_inj and ΔPP_inj versus the same 𝒩_G reference definition as in Section 2.10.

Gate success requires both conditions:

1. CI95(Δθ_inj) lower bound > 0
2. CI95(ΔPP_inj) lower bound > 0

If the gate fails, strong confirmation of Gaussianity is not permitted.

---

# 7. Negative Control (HMDIFF Frozen)

HMDIFF is required.

Discovery cannot be declared if HMDIFF satisfies the discovery criterion defined in Section 3.

Half-mission split products are hard-fixed.

SMICA half missions:

- COM_CMB_IQU-smica_2048_R3.00_hm1.fits
- COM_CMB_IQU-smica_2048_R3.00_hm2.fits

NILC half missions:

- COM_CMB_IQU-nilc_2048_R3.00_hm1.fits
- COM_CMB_IQU-nilc_2048_R3.00_hm2.fits

Local path convention:

- data/external/planck_pr3/<file_name>.fits

SHA256 hashes must be recorded in committed ledgers:

- preregistration/external_inputs/sha256_hm_smica_hm1.txt
- preregistration/external_inputs/sha256_hm_smica_hm2.txt
- preregistration/external_inputs/sha256_hm_nilc_hm1.txt
- preregistration/external_inputs/sha256_hm_nilc_hm2.txt

Binding load call (HM splits):

- HM1 and HM2 are loaded as:
  hp.read_map(path, field=0, nest=False, dtype=numpy.float64, verbose=False)

HMDIFF construction (binding):

- HMDIFF = 0.5 * (HM1 - HM2)

Then apply the same binding order:

1. Smooth HMDIFF
2. Mask to 0
3. Patch
4. Encode
5. Aggregate and compute effects

---

# 8. Environment lock (binding)

Exact numerical reproduction is defined with respect to the frozen environment lock file:

- preregistration/external_inputs/requirements_lock_wsl_ubuntu.txt

SHA256 must be recorded in:

- preregistration/external_inputs/sha256_env_lock_wsl_ubuntu.txt

Primary manuscript values must be generated under an environment that matches this lock.

Runs under other environments may be executed for engineering checks, but their numbers are not admissible for manuscript reporting.

---

# 9. Artifact Policy

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
- all seeds and derived per-realization seeds
- all external input file paths and SHA256 hashes
- patch center list SHA256, manifest SHA256, generator SHA256

---

# 10. Scope Discipline

Restricted to temperature maps.

Excluded:

- polarization
- cosmological parameter inference
- bispectrum fitting
- ontological interpretation without statistical separation

---

# 11. Final Statement

Version 3.1:

- Hard-fixes the two data products SMICA and NILC by name, path convention, and SHA256 ledgers.
- Makes data smoothing a binding pipeline step and fixes its order relative to masking and patching.
- Freezes baseline_seed, null_seed, and inj_seed as explicit integers.
- Freezes 𝒩_G synthesis by specifying the exact healpy.synfast call and the RNG seeding rule.
- Freezes quantile computation by fixing the NumPy quantile method.
- Freezes decision rules and the multi-product conjunction rule.
- Freezes patch array orientation, byte-stream mapping, and block-mean segmentation.
- Freezes injection normalization order and injection seeding.
- Makes the environment lock binding for admissible manuscript numbers.

This document is pre-registration sealed.