---
tags:
  - statistics
  - definition
  - lensfluxanomaly
---
# Statistic Definitions

See [[01_Physics_Definition]] for the physics motivation.

## R_min
```
R_min = min_{|r_i - r_j| < δr}  |F_i - F_j| / (F_i + F_j)
```
- δr = 0.2 arcsec
- r_i = distance from lens centroid (geometric mean of image positions)
- No lens model required
- R_min ∈ [0, 1]

## R_fold
```
R_fold = |F_a - F_b| / (F_a + F_b)
```
- a, b: the two images with positive parity (μ > 0)
- Requires parity from lens model
- R_fold ∈ [0, 1]

## R_cusp
```
R_cusp = |Σ p_i · F_i| / Σ |F_i|
```
- i iterates over the three closest images (minimum sum of pairwise distances)
- p_i = parity of image i (±1)
- For unnormalized simulated data: can use signed magnifications directly
- R_cusp ∈ [0, 1]

## Implementation
- `src/statistic.py`: `compute_rfold()`, `compute_cusp_statistic()`
- `src/compute_rmin.py`: `compute_rmin()`
