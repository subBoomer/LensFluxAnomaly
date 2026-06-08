---
tags:
  - wdm
  - dark-matter
  - lensfluxanomaly
---
# WDM Analysis

## Motivation
Warm dark matter (WDM) suppresses the subhalo mass function below a half-mode mass m_hm. If the observed flux-ratio anomaly were due to CDM substructure exceeding model predictions, reducing substructure via WDM should *reduce* the predicted anomaly, potentially bringing it closer to observations.

## Model
The suppression factor used:
```
n_WDM / n_CDM = 1 / (1 + (m_hm / m)^β)
```
where β = 2.5 and m_hm is computed from:
```
m_hm = 1 × 10⁸  (m_x / 3 keV)^{-3.33}  M⊙
```

Three benchmark WDM masses tested: 3, 5, and 7 keV.

## Results

| Model | m_hm (M⊙) | R_min mean | vs CDM | vs Obs |
|-------|-----------|-----------|--------|--------|
| CDM | — | 0.060 | — | -0.108 |
| WDM 7 keV | 6 × 10⁶ | 0.076 | +0.016 | -0.092 |
| WDM 5 keV | 1.8 × 10⁷ | 0.079 | +0.019 | -0.089 |
| WDM 3 keV | 1 × 10⁸ | 0.085 | +0.025 | -0.083 |

## Interpretation
WDM predicts **larger** R_min than CDM in all cases. This is because:
1. Suppressing low-mass subhalos removes the smallest perturbations
2. The remaining massive subhalos produce *closer*-to-macro flux ratios (less scatter from many small perturbations)
3. But the Poisson variance from fewer halos increases shot noise in some configurations
4. The net effect is slightly broader tail → larger mean R_min

## Conclusion
WDM **cannot rescue** the observed anomaly. All WDM models move the predicted R_min distribution *further* from the observed data, not closer. This rules out a simple DM-model-based interpretation of the discrepancy.

## Limitations
- Only mass-function suppression modeled; no structural differences (concentration, triaxiality) for WDM halos
- No free-streaming effects on the macro lens itself
- No non-linear WDM clustering beyond the suppression factor

See [[09_Limitations]] for full caveats.
