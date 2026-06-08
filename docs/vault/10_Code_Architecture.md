---
tags:
  - code
  - architecture
  - lensfluxanomaly
---
# Code Architecture

## Project Structure
```
LensFluxAnomaly/
├── config.yaml              # Central configuration
├── data/
│   ├── radio_quads.py       # 8 active radio quads + PG1115 placeholder
│   └── curated_quads.py     # 13 CASTES optical quads
├── src/
│   ├── statistic.py         # R_fold, R_cusp, compare functions
│   ├── compute_rmin.py      # R_min statistic
│   ├── perturbation_model.py  # SimplePerturbationModel (CDM)
│   ├── wdm_model.py         # WDMPerturbationModel
│   ├── selection.py         # passes_selection() detection model
│   ├── noise_model.py       # RadioNoise (Gaussian SNR)
│   ├── population.py        # LensPopulation sampler
│   ├── lens_model.py        # MacroLens (SIE+SHEAR+TNFW+LOS)
│   ├── substructure.py      # SubhaloPopulation, LOSPopulation
│   ├── comparison.py        # KS, AD, Wasserstein tests
│   ├── inference.py         # KDE emulator + f_sub posterior
│   ├── config.py            # Config validation
│   └── catalog_utils.py     # build_unified_catalog()
├── run_*.py                 # Entry-point scripts
├── outputs/                 # npz results, png plots
└── docs/
    ├── CONCLUSIONS.md
    ├── EXECUTION_PLAN.md
    ├── PROJECT_GOAL.md
    └── vault/               # Obsidian vault
```

## Run Order
1. `run_rmin_analysis.py` — core R_min analysis (fast)
2. `run_rfold_rcusp_analysis.py` — R_fold/R_cusp analysis (fast)
3. `run_phase_f.py` — WDM comparison (~1 min)
4. `run_validate_model.py` — lenstronomy cross-validation (slow, ~hours)
5. `run_jackknife.py` — jackknife for R_min
6. `run_jackknife_rfrc.py` — jackknife for R_fold/R_cusp
7. `run_outlier_analysis.py` — per-system outlier diagnostics

## Key Results Files
- `outputs/rfold_rcusp_analysis.npz` — R_fold, R_cusp, observed + CDM
- `outputs/phase_b_validation.npz` — lenstronomy cross-validation

## Dependencies
- numpy, scipy, matplotlib
- astropy (cosmology)
- lenstronomy (ray tracing pipeline)
- colossus (halo concentration)
- pyyaml (configuration)
