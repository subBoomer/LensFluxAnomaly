---
tags:
  - model
  - lambda-cdm
  - lensfluxanomaly
---
# ΛCDM Model

Two forward models are used to generate the ΛCDM null distribution:

- [[04_Forward_Model]] — SimplePerturbationModel (phenomenological kernel, 50k realizations)
- [[05_Lenstronomy_Pipeline]] — Full SIE+TNFW+LOS ray tracing via lenstronomy (743 realizations)

## Parameter Values
All fixed at physically motivated defaults. No free parameters fitted to data. Parameters from `config.yaml` and hardcoded defaults in `src/`.

## CDM Subhalo Model
See `src/substructure.py`:
- dN/dm ∝ m^{-1.9}
- Mass range: [10^6, 10^{10}] M⊙
- Spatial: uniform in radius, uniform in angle
- Concentration: Duffy08 (via Colossus)

## WDM Suppression
See [[08_WDM_Analysis]].
