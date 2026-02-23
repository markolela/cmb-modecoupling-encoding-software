# Preregistration Directory

This directory contains the frozen, pre-registered study protocol for the project:

**Beyond Gaussian Sufficiency:  
A Pre-Registered Encoding Test of Scale-Dependent Mode Coupling in Planck 2018 Temperature Maps**

---

## Contents

- `CMB_ModeCoupling_ProjectPlan_v2.9_sealed.md`  
  Final sealed pre-registration document (version 2.9).

- `SHA256.txt`  
  SHA256 hash of the sealed project plan file for integrity verification.

---

## Status

The file `CMB_ModeCoupling_ProjectPlan_v2.9_sealed.md` defines:

- The complete statistical protocol.
- All operational definitions.
- All decision rules.
- Null model construction.
- Sensitivity validation requirements.
- Attribution layer (FFP ensemble).
- Artifact logging policy.

No modifications to this file are permitted after sealing.

Any future methodological changes must appear in a new versioned file (e.g., v3.0) and must not overwrite the sealed document.

---

## Integrity Verification

To verify integrity:

```bash
sha256sum CMB_ModeCoupling_ProjectPlan_v2.9_sealed.md
````

Compare the output with the value stored in `SHA256.txt`.

If the hash differs, the document has been modified.

---

## Relation to Code and Results

* The implementation in `/scripts` must strictly comply with the sealed protocol.
* All results in `/results` must map directly to definitions in the sealed plan.
* Any deviation must be documented explicitly in a new versioned plan.

---

## Scope Discipline

This preregistration applies exclusively to:

* Planck 2018 SMICA and NILC temperature maps.
* The fixed scale set {1,2,4}.
* The primary Gaussian null and fixed FFP10 attribution layer.

No extensions (e.g., polarization, alternative null families, additional scales) fall under this preregistration unless documented in a new sealed version.
