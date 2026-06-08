# Lambda-CDM Model

## Overview
We generate a null distribution of R_min values from a minimal Lambda-CDM forward model. The model simulates perturbations from subhalos using a phenomenological kernel rather than full ray tracing.

## Subhalo Population
- Mass function: dN/dm ∝ m^{-1.9}
- Mass range: 10^6 - 10^9 Msun
- Spatial distribution: exponential radial profile (scale ~1.5 arcsec) centered on lens
- Mean number per system: 200 (drawn from Poisson distribution)

## Perturbation Kernel
For each image i:

F_i_obs = F_i_0 * (1 + sum_k delta_ik)

where:
- F_i_0 = 0.25 (equal fluxes for smooth symmetric quad)
- delta_ik = epsilon * m_k / (r_ik^2 + r_c^2)
- epsilon = 3e-9 (global scaling, tuned to keep perturbations small)
- r_c = 1 kpc (~0.25 arcsec at typical lens redshift)
- r_ik = projected distance from image i to subhalo k

## Noise Model
F_obs = F_perturbed * (1 + N(0, sigma^2))
- sigma = 0.05 (5% flux noise, absorbing microlensing, dust, photometric error)

## Generated Quad Positions
Abstracted SIE-like configurations:
- Einstein radius: lognormal distribution, mean ~1 arcsec
- Ellipticity: uniform [0.1, 0.5]
- Random position angle and orientation
- 4 images placed at radii ~ theta_E * (1 +/- ellipticity_factor)

## Assumptions
1. Subhalos dominate the perturbation signal (no explicit LOS component)
2. The perturbation kernel epsilon * m / (r^2 + r_c^2) captures the leading-order effect
3. Equal macro fluxes (F_i_0 = 0.25) — no intrinsic flux asymmetry
4. Subhalo spatial distribution follows the dark matter halo profile
5. 5% Gaussian noise captures all observational contaminants
