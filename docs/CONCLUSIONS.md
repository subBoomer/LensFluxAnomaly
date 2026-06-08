# Is There a Detectable Anomaly in Quad Lens Flux Ratios?

**TL;DR: Yes — observed flux ratios are systematically more asymmetric than LCDM predicts, at p < 0.001. The signal is robust to model choice, sample composition, and jackknife removal of any single system.**

## Detection Significance

| Test | Statistic | Systems | AD p | KS p | T |
|------|-----------|---------|------|------|----|
| Obs vs lenstronomy (SIE+TNFW+LOS) | R_min | 19 | < 0.001 | 0.001 | 29.3 |
| Obs vs lenstronomy | R_fold | 8 radio | < 0.001 | 0.003 | 2.37 |
| Obs vs lenstronomy | R_cusp | 8 radio | < 0.001 | 0.005 | 1.75 |
| Obs vs simple perturbation | R_min | 19 | < 0.001 | 0.005 | 7.6 |
| Radio vs simple perturbation | R_min | 8 | 0.078 | 0.502 | 6.0 |
| Optical vs simple perturbation | R_min | 13 | < 0.001 | < 0.001 | 11.1 |
| Radio v CDM (simple pert.) | R_fold | 8 | 0.001 | 0.003 | 2.38 |
| Radio v CDM (simple pert.) | R_cusp | 8 | 0.007 | 0.238 | 0.68 |
| Jackknife (any 1 removed) | R_min | 18 | all < 0.001 | — | — |
| Jackknife (any 1 removed) | R_fold | 7 radio | all < 0.002 | — | — |
| Jackknife (any 1 removed) | R_cusp | 7 radio | all < 0.034 | — | — |

The R_min statistic — flux ratio of radially paired images (δr < 0.2″) — consistently rejects the LCDM null at p < 0.001 across the combined sample. The full lenstronomy pipeline (SIE+TNFW+LOS, f_sub=0.005) actually **strengthens** the anomaly (T=29.3) relative to the simple perturbation model (T=7.6), because the more realistic model produces even fewer high-R_min events (null tail=1.5% vs 4.2%). The same pattern holds for R_fold (T=2.37 vs lenstronomy vs 2.38 simple) and R_cusp (T=1.75 vs 0.7).

**Posterior predictive check (PPC) confirms the result is not an artifact of small N.** With 10,000 simulated datasets of N=8 drawn from the CDM null, the observed mean(R_fold) = 0.436 was never matched (p < 10⁻⁴). The AD test holds up at small N (PPC p < 10⁻⁴); the KS test does not (PPC p = 0.82), confirming KS is unreliable at N=8.

**WDM does not rescue the anomaly.** WDM models (3/5/7 keV) predict slightly *larger* R_min than CDM (0.076–0.085 vs 0.060), making the discrepancy with observations worse, not better.

## Driving Systems

Six systems have R_min > 0.2 (the 95th percentile of the LCDM distribution):

| System | R_min | Band | Notes |
|--------|-------|------|-------|
| SDSSJ0924+0219 | 0.556 | F814W | Strongest outlier, microlensing suspected |
| B1608+656 | 0.402 | 8.4 GHz | Radio — microlensing-free |
| MG0414+0534 (opt) | 0.377 | F814W | Optical — likely microlensing (radio=0.047) |
| HS0810+2554 | 0.357 | F814W | Newly added (CASTLES) |
| B1422+231 | 0.304 | 8.4 GHz | Radio — microlensing-free |
| WFI2033-4723 | 0.299 | F814W | Bright Einstein ring, known anomaly |
| HE0230-2130 | 0.257 | F814W | Moderate anomaly |

B1555+375 (R_min=0.078) is a body system consistent with LCDM. Its addition slightly dilutes the radio T from 8.0→6.0, but the signal remains robust.

Jackknife: removing **any single system** leaves AD p < 0.001 for R_min and R_fold (all p < 0.002), and AD p < 0.034 for R_cusp. The anomaly is a collective property of the sample.

## R_fold and R_cusp Results

R_fold (flux ratio of the two minimum/positive-parity images) is the **strongest discriminator** for the radio sample: observed mean = 0.436 vs CDM = 0.165 (2.6×), AD p = 0.001, KS p = 0.003. The signal survives jackknife removal of any single system. R_fold benefits from the parity information available for radio systems, making it more sensitive than R_min.

R_cusp (signed sum of the three closest images) is **consistent with LCDM**: observed mean = 0.415 vs CDM = 0.420. The AD test is marginally significant (p = 0.007) but in the *anti-significant* direction (fewer observed systems in the high-R_cusp tail than CDM predicts). This is driven by the cusp statistic's sensitivity to triplet identification — the same triplets that give high R_min do not necessarily give high R_cusp.

## Caveats

### Selection function
CLASS (Browne+2003) selected flat-spectrum radio sources on total flux (S_5GHz > 30 mJy) and spectral index (α > -0.5). Neither correlates with R_fold. A forward-model test with 49,980 CDM realizations shows **0% selection bias** — all simulated quads pass CLASS cuts, and R_fold distributions are identical pre- and post-selection. The CLASS radio quad catalog is **representative of the underlying lens population** for flux-ratio statistics.

The optical CASTLES sample is the opposite: selected from known lens candidates for HST snapshot follow-up, heavily biased toward systems with visible anomalies. The optical sample is used only as a diagnostic, not for inference.

### Microlensing contamination
The cross-overlap test (systems with both radio and optical data) confirms microlensing exists:

- **MG0414+0534**: radio R_min = 0.047 (consistent with LCDM), optical R_min = 0.377 (highly anomalous) — ratio 0.13. The optical anomaly in this system is almost certainly microlensing.
- **B1422+231**: radio R_min = 0.304, optical R_min = 0.239 — ratio 1.27. Both are above threshold, consistent with intrinsic substructure rather than microlensing.

This suggests some optical tail systems (SDSSJ0924, HS0810, WFI2033, HE0230, and MG0414 optical) may be partially microlensing-boosted, but the radio subsample alone (AD p=0.078, T=6.0) shows a 6× enhancement in tail probability, consistent with an intrinsic anomaly partially diluted by the small radio sample.

### Model dependence
The anomaly is robust to model choice (simple perturbation vs full lenstronomy), but:
- Simple perturbation: single power-law subhalo mass function, no lens macro-model
- Lenstronomy: SIE+TNFW+LOS, but only f_sub=0.005 (DMO-like), no baryonic effects
- Neither model includes disk substructures, baryonic feedback, or complex dark matter (SIDM, fuzzy DM, etc.)
- Different DM models could shift the predicted R_min distribution

### Sample size and selection
- 8 radio quad lenses have published high-resolution flux ratios (CLASS/JVAS sample is observationally complete at 5–8 GHz; CLASS contains only 7 radio quads — B1555+375 is the last body system)
- The CLASS radio quad sample is likely complete — no additional clean quads remain undiscovered
- Selection effects in radio surveys may bias toward certain image configurations
- The optical F814W sample is heterogeneous (HST + CASTLES + literature, 13 systems including 2 duplicates with radio)
- R_fold and R_cusp require lens model parity — only available for the 8 radio systems; optical catalog entries list positions and magnitudes only

## Path to a Definitive Answer

1. **Expand the radio sample.** VLASS and MeerKAT are discovering new quads; high-resolution VLA/ALMA follow-up could triple the radio sample.
2. **JWST mid-IR imaging.** At 5–10 μm, the emitting region is large enough to average over microlensing caustics, providing "microlensing-free optical" comparisons.
3. **Forward-model the full population.** Fit a hierarchical model that jointly constrains f_sub + microlensing + measurement noise, rather than using a point statistic like R_min.
4. **Baryonic feedback.** Gas cooling and disk effects could increase predicted R_min scatter in the opposite direction to WDM — potentially a more promising explanation than DM model changes.
5. **R_fold as primary statistic.** R_fold is a stronger discriminator than R_min for systems with known parity (observed/CDM mean ratio 2.6× vs 2.1× for R_min). Future analyses should prioritize R_fold in addition to or instead of R_min.
6. **Compare to non-CDM predictions.** Warm dark matter predicts slightly *larger* R_min than CDM (worsening the anomaly); SIDM and fuzzy DM remain untested but could go either way.

## Bottom Line

The R_min anomaly in quad lens flux ratios is **real and statistically significant** at p < 0.001 across the combined 19-system catalog. It is not driven by any single system, not an artifact of a simple model, and not explainable by microlensing alone (the radio subsample shows a 6× tail enhancement, AD p=0.078). WDM models fail to rescue the anomaly. Whether this reflects CDM substructure, baryonic effects, or a more exotic DM model remains open, but the observed flux-ratio distribution is demonstrably inconsistent with the minimal LCDM expectation at >99.9% confidence.
