---
tags:
  - lenstronomy
  - model
  - lensfluxanomaly
---
# Lenstronomy Pipeline

## Purpose
Cross-validate the simple perturbation model against full SIE+TNFW+LOS ray tracing. If the simple model underestimates the scatter, the more realistic pipeline will produce *fewer* extreme events, making the observed anomaly *stronger*.

## Components

### Macro Lens (SIE + External Shear)
- Singular Isothermal Ellipsoid (SIE)
- External shear: gamma_ext drawn from exponential(scale=0.05)
- Axis ratio q: truncated normal (mean=0.8, sigma=0.15), clipped [0.3, 1.0]
- Velocity dispersion sigma_v: log-normal (mu=ln(250), sigma=0.15)

### Lens Population
- Lens redshift z_l: comoving-volume-weighted dV/dz * (1+z)²
- Source redshift z_s: uniform over [1.0, 3.0]
- Source position: sampled within 0.8× caustic to ensure 4 images

### Line of Sight (LOS)
- Correlated convergence + shear from a multivariate normal
- sigma_kappa = 0.030 × sqrt(D_C(z_s) / D_C(2)) × (1+z_l)/1.5
- sigma_gamma = 0.8 × sigma_kappa, corr(kappa, gamma1) = 0.4

### Subhalo Population (TNFW)
- Same dN/dm ∝ m^{-1.9} as simple model, over [10^6, 10^{10}] M⊙
- Truncated NFW profile via lenstronomy
- Concentration from Colossus Duffy08
- f_sub = 0.005 (DMO equivalent)

### Selection
- Flux limit: total lensed flux ≥ 1 mJy
- Resolution: sigmoid probability below 0.25 arcsec
- Quad detection probability: 0.7

### Noise
- Gaussian noise at SNR = 50

## Output
743 quad realizations from 2000 seeds (37% selection efficiency) at f_sub = 0.005.

## Key Results
| Statistic | Obs mean | Lenstronomy mean | KS p | AD sig | Tail T |
|-----------|---------|-----------------|------|--------|--------|
| R_min | 0.168 | 0.047 | 0.001 | < 0.001 | 29.3 |
| R_fold | 0.436 | 0.164 | 0.003 | < 0.001 | 2.37 |
| R_cusp | 0.415 | 0.162 | 0.005 | < 0.001 | 1.75 |

The lenstronomy null is tighter than the simple model null in all three statistics, confirming that the more realistic model **strengthens** the anomaly.
