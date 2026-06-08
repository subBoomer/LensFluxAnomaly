# Lens Flux-Ratio Anomaly Statistics (Minimal-Model ΛCDM Test)

## Objective
Build a reproducible computational framework that tests whether observed flux-ratio asymmetries in
Hubble quad gravitational lenses deviate from a minimal ΛCDM-consistent expectation.

**Statistic**: R_min = min_{|r_i - r_j| < δr} |F_i - F_j| / (F_i + F_j)

**Process**:
1. Compute empirical R_min distribution from observed quads
2. Generate ΛCDM mock distribution using a simplified subhalo perturbation model
3. Compare distributions statistically (KS, Anderson-Darling, tail ratio)
4. Identify presence or absence of heavy-tail excess

**Core constraint**: No per-system lens fitting. No qualitative classification.

## Original Prompt (master prompt)

CODING AGENT MASTER PROMPT
Project: Lens Flux-Ratio Anomaly Statistics (Minimal-Model ΛCDM Test)

### 1. Objective

Build a reproducible computational framework that tests whether observed flux-ratio asymmetries in Hubble quad gravitational lenses deviate from a minimal ΛCDM-consistent expectation.

This is done using a model-light observable statistic:

R_min = min_{ (i,j): |r_i - r_j| < δr_min } |F_i - F_j| / (F_i + F_j)

Goal:
- Compute empirical distribution of R_min
- Generate ΛCDM mock distribution using a simplified subhalo perturbation model
- Compare distributions statistically
- Identify presence or absence of a heavy-tail excess

No per-system lens fitting is allowed.

### 2. Data Sources (public only)

Primary lens catalogs — use ALL of the following:
- HST strong lens systems
- Hubble Legacy Archive (HLA)
- CASTLES lens sample
- SLACS survey (Sloan Lens ACS Survey catalog: positions + fluxes)
- BELLS survey (BOSS Emission-Line Lens Survey)
- DES / KiDS strong lens catalogs
- Quad lens subsets only

### 3. Data extraction rules

For each system:
- Extract: image positions (θ_x, θ_y), fluxes F_i in a single consistent band (or nearest filter match), lens centroid θ_lens
- Normalize: F_i → F_i / ΣF_i
- Filter: ONLY quad lenses; discard systems with missing fluxes or ambiguous image IDs

### 4. Core observable computation

Step 1 — Radial coordinate: r_i = ‖θ_i − θ_lens‖
Step 2 — Pairing rule (no lens model): pairs with |r_i − r_j| < δr, default δr = 0.2 arcsec
Step 3 — Statistic: R_min = min |F_i - F_j| / (F_i + F_j)

Store: per-system R_min, full population distribution.

### 5. ΛCDM null model (minimal forward simulation)

#### 5.1 Macro lens baseline
Assume elliptical isothermal lens (SIE-like behavior, but NOT explicitly modeled) producing symmetric quad image positions (abstracted). Do NOT fit real lenses.

#### 5.2 Subhalo population
For each lens: draw subhalos with dN/dm ∝ m^−1.9, mass range 10^6−10^9 M⊙, spatial distribution NFW-tracing host halo.

#### 5.3 Flux perturbation model (critical simplification)
For each image: F_i^obs = F_i^0 (1 + Σ_k δ_ik)
where δ_ik = ε · m_k / (r_ik² + r_c²)

Parameters: ε (global scaling, fit to keep perturbations small), r_c (softening scale, ~1 kpc equivalent).

This is a proxy perturbation kernel, NOT full ray tracing.

#### 5.4 Output
Generate simulated R_min distribution, multiple realizations (N ≥ 10,000 systems).

### 6. Contaminant handling (simplified)

Flux noise: F_i → F_i (1 + N(0, σ²)), default σ = 0.05.
This absorbs microlensing, dust, photometric error.

### 7. Statistical comparison

- Empirical CDF vs simulated CDF
- KS test
- **Anderson-Darling test (preferred)**
- Tail ratio: T = P_obs(R_min > 0.2) / P_ΛCDM(R_min > 0.2)

Primary signal: excess in high-R_min tail.

### 8. Software stack

Python 3.11+, numpy, scipy, pandas, astropy, matplotlib, scikit-learn (optional),
astroquery (for catalog ingestion). No proprietary tools.

### 9. Project structure

```
LensFluxAnomaly/
├── 00_Project_Overview.md
├── 01_Physics_Definition.md
├── 02_Data_Sources.md
├── 03_Statistic_Definition.md
├── 04_LambdaCDM_Model.md
├── 05_Data_Pipeline.md
├── 06_Simulation_Model.md
├── 07_Analysis_Methods.md
├── 08_Results.md
├── 09_Limitations.md
├── 10_Code_Architecture.md
├── data/
│   ├── raw/ | processed/ | catalogs/
├── src/
│   ├── ingest.py | compute_statistic.py | simulation.py | lens_filters.py | stats.py
├── notebooks/
│   ├── exploratory_analysis.ipynb | simulation_tests.ipynb
├── outputs/
│   ├── distributions/ | plots/
└── README.md
```

### 10. Required Obsidian note contents

- 00_Project_Overview.md: objective, hypothesis, constraints (no lens refitting rule)
- 03_Statistic_Definition.md: full derivation of R_min, model-independence explanation
- 06_LambdaCDM_Model.md: full Monte Carlo description, assumptions explicitly listed
- 08_Results.md: only statistical outputs, no interpretation until comparison stage

### 11. Execution phases

1. Data ingestion — collect catalogs, unify format, extract quad systems
2. Observable computation — compute R_min, store dataset
3. ΛCDM simulation — generate mock lens population, compute simulated distribution
4. Statistical comparison — KS test, tail analysis, distribution plots
5. Anomaly evaluation — quantify deviation strength, identify outlier systems

### 12. Final expected output

- Empirical R_min distribution
- ΛCDM simulated distribution
- Statistical comparison metrics
- List of strongest anomalous systems (if any)
- Reproducible codebase
- Fully structured Obsidian research vault

### 13. Non-negotiable constraints

- No per-system lens modeling
- No qualitative classification
- No alternative statistics
- No external API dependence after ingestion stage
- Everything must be reproducible offline once data is downloaded
