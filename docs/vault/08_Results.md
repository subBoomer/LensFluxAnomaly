---
tags:
  - results
  - lensfluxanomaly
---
# Results

All results from the final 19-system catalog (8 radio + 13 optical, 2 duplicates deduplicated). CDM null from 50,000 SimplePerturbationModel realizations. Lenstronomy null from 743 SIE+TNFW+LOS realizations.

## R_min — Combined Sample (19 systems)

| Metric | Value |
|--------|-------|
| Observed mean | 0.168 ± 0.152 |
| CDM simple model mean | 0.060 ± 0.064 |
| KS p | 0.005 |
| AD sig | < 0.001 |
| Tail ratio T | 7.6 |
| Observed tail (R > 0.2) | 6/19 (32%) |
| Simulated tail | 4.2% |

## R_min — Radio Sample (8 systems)

| Metric | Value |
|--------|-------|
| Observed mean | 0.128 ± 0.135 |
| CDM simple model mean | 0.060 ± 0.064 |
| KS p | 0.502 |
| AD sig | 0.078 |
| Tail ratio T | 6.0 |
| Tail (R > 0.2) | 2/8 (B1608, B1422) |

## R_min — Optical Sample (13 systems)

| Metric | Value |
|--------|-------|
| Observed mean | 0.214 ± 0.152 |
| CDM simple model mean | 0.060 ± 0.064 |
| KS p | < 0.001 |
| AD sig | < 0.001 |
| Tail ratio T | 11.1 |
| Tail (R > 0.2) | 6/13 (SDSSJ0924, MG0414_opt, HS0810, WFI2033, B1422_opt, HE0230) |

## R_fold — Radio Sample (8 systems, primary statistic)

| Metric | Value |
|--------|-------|
| Observed mean | 0.436 ± 0.260 |
| CDM simple model mean | 0.165 ± 0.138 |
| KS p | 0.003 |
| AD p (permutation) | 0.001 |
| Tail ratio T | 2.38 |

## R_cusp — Radio Sample (8 systems)

| Metric | Value |
|--------|-------|
| Observed mean | 0.415 ± 0.300 |
| CDM simple model mean | 0.420 ± 0.144 |
| KS p | 0.238 |
| AD p (permutation) | 0.007 |
| Tail ratio T | 0.68 |

R_cusp is consistent with or anti-significant relative to ΛCDM.

## Lenstronomy Cross-Validation (743 realizations, f_sub=0.005)

| Statistic | Obs mean | Lenstronomy mean | KS p | AD sig | Tail T |
|-----------|---------|-----------------|------|--------|--------|
| R_min (N=19) | 0.168 | 0.047 | 0.001 | < 0.001 | 29.3 |
| R_fold (N=8) | 0.436 | 0.164 | 0.003 | < 0.001 | 2.37 |
| R_cusp (N=8) | 0.415 | 0.162 | 0.005 | < 0.001 | 1.75 |

Lenstronomy produces fewer extreme events than the simple model (tighter null) → anomaly is **stronger** against the more realistic model.

## Jackknife

Removing any single system leaves AD p < 0.001 for R_min (18 systems) and AD p < 0.002 for R_fold (7 radio). The anomaly is a collective property of the sample, not driven by any single outlier.

## WDM Analysis

| Model | Half-mode mass | R_min mean | R_min vs CDM |
|-------|---------------|-----------|-------------|
| CDM | — | 0.060 | — |
| WDM 7 keV | 6 × 10⁶ M⊙ | 0.076 | Larger (worse) |
| WDM 5 keV | 1.8 × 10⁷ M⊙ | 0.079 | Larger (worse) |
| WDM 3 keV | 1 × 10⁸ M⊙ | 0.085 | Larger (worse) |

WDM cannot rescue the anomaly — all models predict *larger* R_min than CDM, moving further from the data.

## Tail Systems (R_min > 0.2)

| System | R_min | Band | Notes |
|--------|-------|------|-------|
| SDSSJ0924+0219 | 0.556 | F814W | Strongest — microlensing suspected |
| B1608+656 | 0.402 | 8.4 GHz | Radio — microlensing-free |
| MG0414+0534 | 0.377 | F814W | Optical — radio version is 0.047 (microlensing) |
| HS0810+2554 | 0.357 | F814W | CASTLES addition |
| B1422+231 | 0.304 | 8.4 GHz | Radio — microlensing-free |
| WFI2033-4723 | 0.299 | F814W | Known Einstein ring system |
| HE0230-2130 | 0.257 | F814W | Moderate anomaly |
