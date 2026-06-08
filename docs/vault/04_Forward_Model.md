---
tags:
  - model
  - simulation
  - lensfluxanomaly
---
# Forward Model — Simple Perturbation

## Philosophy
A minimal model for what ΛCDM predicts for flux ratios, without ray tracing. Used to generate the null distribution for R_min, R_fold, and R_cusp.

## Abstract Quad Geometry
Each mock system generates 4 lensed images with:
- Einstein-ring-like positions: angles [0, π/2, π, 3π/2], radii [r1, r2, r1, r2]
- r1, r2 drawn from distributions matching observed CLASS quad separations
- Parity assignment: [1, -1, 1, -1] (alternating around the ring)
- Equal macro fluxes: F_i^0 = 0.25 (normalized)

## Subhalo Population
- Mass function: dN/dm ∝ m^{-1.9} over [10^6, 10^9] M⊙
- Spatial: uniform in radius within θ_E, uniform in angle
- Mean number: n_sub = 200 (Poisson drawn per realization)
- No tidal stripping or spatial bias modeled

## Perturbation Kernel
Each subhalo perturbs each image:
```
δ_ik = ε · m_k / (r_ik² + r_c²)
```
where ε = 5 × 10^{-9} (amplitude), r_c = 1 kpc (softening).

Observed flux: F_i^obs = F_i^0 (1 + Σ_k δ_ik)

## Noise Model
```
F_i^final = F_i^obs × N(0, σ_noise)
```
where σ_noise = 0.05 (5% multiplicative scatter).

## Parameters
| Parameter | Value | Notes |
|-----------|-------|-------|
| ε | 5 × 10^{-9} | Perturbation amplitude |
| r_c | 1 kpc | Softening radius |
| m_min | 10^6 M⊙ | Min subhalo mass |
| m_max | 10^9 M⊙ | Max subhalo mass |
| α | 1.9 | Mass function slope |
| n_sub | 200 | Mean subhalo count |
| σ_noise | 0.05 | Flux scatter |

## WDM Variant
Implemented in `WDMPerturbationModel`: same kernel, but subhalo mass function is suppressed:
```
n_WDM(m) / n_CDM(m) = 1 / (1 + (m_hm/m)^β)
```
where β = 2.5, m_hm = 10^8 × (m_x/3 keV)^{-3.33} M⊙.

## Coverage
50,000 realizations per model variant (CDM, WDM 3/5/7 keV).
