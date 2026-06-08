---
tags:
  - statistics
  - methods
  - lensfluxanomaly
---
# Statistical Methods

## Primary Test: Anderson-Darling
The Anderson-Darling (AD) k-sample test is the primary comparison. It is sensitive to differences in the tails of distributions, which is where the flux-ratio anomaly appears (high-R_min tail).

AD p-values are computed via:
1. Asymptotic approximation (scipy default)
2. Permutation test (500-2000 resamples) for small-sample calibration

## Secondary Test: Kolmogorov-Smirnov
The KS two-sample test is reported as a secondary diagnostic. It is less sensitive to tail differences than AD, so it is expected to give larger p-values.

## Tail Ratio T
```
T = P_obs(R > 0.2) / P_sim(R > 0.2)
```
T is a descriptive metric, not a formal test statistic. It depends on the threshold choice (0.2 is the ~95th percentile of the CDM null distribution).

## Jackknife Leave-One-Out
Each observed system is removed one at a time, and the AD test is recomputed. This tests whether the signal is driven by a single outlier.

## Permutation Test
For the AD test at small sample sizes (N=8), a permutation p-value is computed:
1. Combine observed and simulated arrays
2. Shuffle, split into two groups of same sizes
3. Compute AD statistic
4. Repeat 500-2000 times
5. p = fraction of permutations with AD ≥ observed AD

## Posterior Predictive Check (PPC)
Planned but not yet implemented. Would replace AD p-values with calibrated significance:
1. Simulate N datasets at best-fit model parameters
2. Compute test statistic T for each
3. Compare observed T to the PPC distribution of T
4. Calibrated p-value = fraction of PPC draws with T ≥ observed T

## Current Limitations
- No likelihood function defined → no Bayesian inference on f_sub
- No joint test over multiple statistics → risk of p-value inflation
- No selection function calibration → risk of selection-dominated inference

See [[09_Limitations]] for full details.
