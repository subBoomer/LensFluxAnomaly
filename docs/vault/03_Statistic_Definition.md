# Statistic Definition

## R_min

R_min = min_{ (i,j): |r_i - r_j| < delta_r } |F_i - F_j| / (F_i + F_j)

where:
- r_i = sqrt((x_i - x_c)^2 + (y_i - y_c)^2) is the radial distance from the lens centroid
- x_c, y_c = geometric center of image positions (mean of x_i, y_i)
- F_i are normalized fluxes (F_i = F_raw_i / sum(F_raw))
- delta_r = 0.2 arcsec (default)

## Rationale
This statistic is designed to be:
- **Model-independent**: requires no lens model, no parity information, no redshifts
- **Sensitive to fold-pair anomalies**: pairs images at similar radii, which approximates the fold pair selection
- **Normalized**: |F_i - F_j|/(F_i + F_j) is bounded between 0 and 1

## Interpretation
- R_min = 0: the closest pair in radius has identical fluxes (smooth lens)
- R_min -> 1: extreme flux asymmetry in the closest radial pair
- High R_min values indicate strong perturbations from substructure

## Comparison Statistics
- **Anderson-Darling test** (primary): sensitive to tail differences between distributions
- **KS test** (secondary): detects general distribution differences
- **Tail ratio T**: P_obs(R_min > 0.2) / P_sim(R_min > 0.2) — quantifies excess in the high-R_min tail
