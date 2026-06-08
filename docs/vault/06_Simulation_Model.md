# Simulation Model

## Overview
The simulation generates synthetic quad lens systems with perturbed flux ratios using a phenomenological subhalo model. This is deliberately simpler than full ray-tracing simulations (no lenstronomy, no macro model fitting).

## Architecture
File: `src/perturbation_model.py`
Class: `SimplePerturbationModel`

## Per-Realization Steps
1. Generate abstract quad positions (SIE-like cross pattern)
2. Draw subhalo population from mass function
3. Compute flux perturbations using kernel
4. Apply Gaussian noise
5. Normalize fluxes
6. Compute R_min

## Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| epsilon | 3e-9 | Perturbation scaling |
| r_c | 1 kpc | Softening scale |
| m_min | 1e6 Msun | Minimum subhalo mass |
| m_max | 1e9 Msun | Maximum subhalo mass |
| alpha | 1.9 | Mass function slope |
| sigma | 0.05 | Flux noise amplitude |
| n_sub_mean | 200 | Mean subhalo count (Poisson) |

## Output
- Array of R_min values from N realizations
- Typical yield: ~99.9% of draws produce valid quads
- Default: 50,000 realizations

## Performance
- ~6,400 realizations/second on single core
- 50k realizations in ~8 seconds
