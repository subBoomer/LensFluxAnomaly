---
tags:
  - project-overview
  - lensfluxanomaly
---
# Project Overview

## Objective
Test whether observed flux-ratio asymmetries in quad gravitational lenses deviate from a minimal ΛCDM-consistent expectation using model-light statistics.

## Hypothesis
ΛCDM predicts a specific level of flux-ratio anomalies from dark matter substructure + line-of-sight structure. If observed anomalies exceed predictions, this could indicate:
- Higher substructure abundance than DMO models predict
- Additional physics beyond minimal ΛCDM
- Systematic biases in the observational data (selection, microlensing)

## Statistics
Three complementary statistics, each with different assumptions and sensitivities:

| Statistic | Target | Data Requirement |
|-----------|--------|-----------------|
| R_min | Radial pairing, model-light | Positions + fluxes, no lens model |
| R_fold | Fold pair (two minima) | Positions + fluxes + parity |
| R_cusp | Three closest images | Positions + fluxes + parity |

## Constraints
- No per-system lens fitting
- Public data only (CLASS, CASTES, literature)
- Reproducible via seeded random number generation

## Current Status
- **19 unique systems** compiled (8 radio + 13 optical, 2 duplicates)
- **R_min**: AD sig < 0.001 across combined sample, T = 7.6
- **R_fold** (primary, 8 radio): AD p = 0.001, 2.6× mean excess
- **R_cusp**: consistent with ΛCDM (AD p = 0.007, anti-significant)
- **Lenstronomy validation**: strengthens all signals
- **WDM**: fails to rescue the anomaly (makes it worse)

See [[08_Results]] for full numbers.
