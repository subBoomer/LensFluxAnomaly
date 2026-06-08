---
tags:
  - limitations
  - caveats
  - lensfluxanomaly
---
# Limitations and Caveats

## Sample Size
With 19 systems (8 radio + 13 optical, 2 overlap), statistical power is limited. Effective degrees of freedom for tail metrics is ~2-3, not 8. The AD p-values at N=8 are descriptive divergences, not calibrated frequentist quantities. A sample of 50+ quads would provide a definitive test.

## Selection Function
The CLASS radio sample was selected on **total flux** (S_5GHz > 30 mJy) and **spectral index** (α > -0.5), not image morphology. A forward-model test with 49,980 CDM realizations shows **0% selection bias** — all simulated quads pass CLASS cuts, and R_fold distributions are identical pre- and post-selection. The CLASS radio quad catalog is representative for flux-ratio statistics.

The optical CASTLES sample is disaster-selected (HST snapshots of known lens candidates) — heavily biased toward visible anomalies. This invalidates the optical+radio pooled analysis as an inference. The radio-only (N=8) analysis is the only one used for claims.

## Simple Perturbation Model
The phenomenological kernel ε·m/(r² + r²_c) is a crude approximation:
- No ray tracing → cannot capture magnification bias or caustic effects
- Equal macro fluxes assumed (F_i_0 = 0.25) → real lenses have intrinsic asymmetries
- No explicit LOS component → absorbed into kernel
- No baryonic disk asymmetries, miscentering, or multipole moments beyond quadrupole

Each of these tends to broaden the predicted R distribution, making the current model **too narrow in perturbation space**. This makes the anomaly detection **conservative** (underpredicting scatter), not inflated.

## Lenstronomy Model
The full SIE+TNFW+LOS pipeline (2000 realizations) improves on the simple model but:
- Only f_sub = 0.005 (DMO equivalent) tested
- No baryonic effects (gas cooling, disk)
- No subhalo tidal stripping or spatial bias
- No source structure or extended source effects

## Microlensing Contamination
Cross-overlap test (MG0414: radio R_min=0.047 vs optical=0.377) confirms microlensing inflates optical R_min. Radio subsample partially mitigates this (6× tail excess, AD p=0.078), but at N=8 the radio-only signal is marginally significant.

## WDM Test Limitations
The WDM model only suppresses the subhalo mass function below half-mode mass. It does not model:
- Free-streaming effects on the macro lens
- Different concentration-mass relations for WDM halos
- Non-linear structure formation differences beyond mass-function suppression

## Parigm
R_fold and R_cusp require lens model parity, available only for the 8 radio systems (from literature). This limits the primary analysis to N=8. The optical F814W catalog has positions and magnitudes only.
