---
tags:
  - simulation
  - model
  - lensfluxanomaly
---
# Simulation Model

See [[04_Forward_Model]] for the simple perturbation model.
See [[05_Lenstronomy_Pipeline]] for the lenstronomy validation.
See [[08_WDM_Analysis]] for WDM predictions.

## Realization Counts
| Model | Realizations | Description |
|-------|-------------|-------------|
| CDM simple | 50,000 | SimplePerturbationModel |
| WDM 3 keV | 50,000 | WDMPerturbationModel |
| WDM 5 keV | 50,000 | WDMPerturbationModel |
| WDM 7 keV | 50,000 | WDMPerturbationModel |
| Lenstronomy | 743 | SIE+TNFW+LOS (3,917 macros with 4 images, 743 pass selection) |

## Random Seeds
All model realizations start from `numpy.random.default_rng(42)`. Each subsequent draw increments the seed. The global random state is never used (all random generation is through local Generator instances).
