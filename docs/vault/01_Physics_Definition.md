---
tags:
  - physics
  - lensfluxanomaly
---
# Physics Definition

## R_min
Flux ratio of radially paired images:
```
R_min = min_{|r_i - r_j| < δr}  |F_i - F_j| / (F_i + F_j)
```
- δr = 0.2 arcsec (standard radial pairing threshold)
- r_i = distance from lens centroid
- No lens model needed — only positions and fluxes

R_min is the most model-independent statistic but has the weakest physical grounding. It is used as a **diagnostic** in this project, not as the primary inference statistic.

## R_fold
Flux ratio of the two minimum (positive-parity) images:
```
R_fold = |F_a - F_b| / (F_a + F_b)
```
- a, b are the two images with μ > 0 (positive magnification)
- Directly tied to fold caustic catastrophe expansion
- Requires parity from lens model (literature for radio systems)
- Most standard in literature (Dalal & Kochanek lineage)

R_fold is the **primary inference statistic** in this project.

## R_cusp
Signed sum ratio of the three closest images:
```
R_cusp = |p_a·F_a + p_b·F_b + p_c·F_c| / (F_a + F_b + F_c)
```
- p_i = parity (±1) of image i
- Triplet is the three images with minimum pairwise distance sum
- Requires parity + positions
- Sensitive to third-order caustic expansion
- More model-dependent than R_fold

R_cusp is analyzed but not used for primary claims.

## Relationship Between Statistics
These three statistics are not independent test statistics under one null model. They differ in:
- Image pairing rule (radial vs parity vs triplet)
- Sensitivity to shear orientation and caustic curvature
- Robustness to mis-centered lens position

A joint test (Mahalanobis distance in (R_fold, R_cusp) space) is possible but not implemented.
