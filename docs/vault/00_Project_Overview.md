# Project Overview

## Objective
Test whether observed flux-ratio asymmetries in Hubble quad gravitational lenses deviate from a minimal Lambda-CDM-consistent expectation.

## Hypothesis
Lambda-CDM predicts a specific level of flux-ratio anomalies from dark matter substructure along the line of sight. If the observed anomalies are stronger than predicted, this could indicate:
- A higher substructure abundance than expected
- Additional physics beyond the standard model
- Systematic biases in the observational data

## Observable Statistic
R_min = min_{|r_i - r_j| < delta_r} |F_i - F_j| / (F_i + F_j)

This is a model-light statistic requiring no per-system lens fitting — only image positions and fluxes relative to the lens centroid.

## Constraints
- No per-system lens modeling
- No qualitative classification
- Reproducible offline computation
- Public data only

## Current Status
- 15 quad lens systems compiled (7 radio + 8 optical)
- R_min: KS p=0.004, AD sig<0.001, tail ratio T=8.0
- Anomaly detected at high significance in the combined sample
- Optical sample drives the signal (AD sig<0.001, T=12.0)
- Radio sample alone borderline (AD sig=0.078, T=6.8)
