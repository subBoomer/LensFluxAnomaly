# LensFluxAnomaly — Project Status & Next Steps

## Project Goal

Constrain the dark matter substructure abundance parameter f_sub (fraction of mass in subhalos ~10⁶–10⁹ M☉) by comparing observed R_fold and R_cusp flux ratio anomalies in quadruple-image gravitational lens systems against a simulated lens population with baryonic + dark matter + substructure components.

The core idea: smooth lens models predict certain ratios of image fluxes. If observed flux ratios deviate systematically, the excess perturbation can be attributed to low-mass dark matter substructure along the line of sight or in the lens halo.

**Key reference methodology:** Dalal & Kochanek (2002), Xu et al. (2015), Nierenberg et al. (2017).

---

## Current Codebase Architecture

```
LensFluxAnomaly/
├── run_inference.py          # Main entry: draws lens systems, computes R_fold/R_cusp, KS test vs observations
├── run_fsub_grid.py          # Runs f_sub grid across values, with --fast parallel mode
├── run_ablation.py           # Runs ablation tests (toggle substructure/LOS/2D vs 3D)
├── run_systematic.py         # Systematic tests (SNR sweep, concentration model sweep)
├── config.yaml               # Central config: n_systems=5000, radio_snr=20, etc.
├── data/
│   └── radio_quads.py        # Observed lens catalog: 8 active systems (7 radio + 1 optical)
├── src/
│   ├── substructure.py       # LOSPopulation (Gaussian κ/γ), SubhaloPopulation (TNFW profiles)
│   ├── inference.py          # KDE emulator + grid integration posterior for f_sub
│   ├── plot_ablation.py      # Ablation comparison plot generator
│   └── run_ablation.py       # Ablation runner logic
└── outputs/                  # All results (npz checkpoints, png plots)
```

---

## What We've Done (Completed Phases)

### Phase C1: LOS Population
- Implemented `LOSPopulation.realise()` in `src/substructure.py` with Gaussian distributions for κ and γ
- σ_κ ≈ 0.02–0.04, calibrated from Hilbert+2009 ray-tracing through N-body simulations
- Scaled by source-lens geometry (comoving distances)
- Integrated into both parallel and sequential paths in `run_inference.py`
- Added LOS to "full" config in `run_ablation.py`

### Phase C2: f_sub Posterior Inference
- 6-value f_sub grid: [0.000, 0.003, 0.005, 0.008, 0.010, 0.020]
- Created `src/inference.py`: KDE-based R_fold distribution emulator + P(f_sub|R_obs) via grid integration
- `run_fsub_grid.py` with --fast parallel mode (14× speedup over serial)
- **Result: f_sub posterior is flat —** the LOS alone already explains the anomaly. Subhalos on top of LOS don't add distinguishable signal.

### Phase C3: Cusp Statistic
- Wired `compute_cusp_statistic` into `run_inference.py`
- Stores both R_fold and R_cusp in checkpoint
- Fixed signed-magnification bug (was using signed fluxes with noise; now uses noise-free signed magnifications for simulated data's cusp formula)

### Task 2: Ablation Plot
- `src/plot_ablation.py` → `outputs/ablation_comparison.png`
- 6 configs × 2 statistics (R_fold, R_cusp) side-by-side with KS p-values
- Configs: smooth, smooth+LOS, smooth+sub, smooth+sub+LOS, full2D, full3D

### Task 3: Systematic Tests
- SNR sweep (20, 50, 100) — results robust across SNR
- Concentration model sweep (duffy08, ishiyama21, diemer15) — results robust
- All show KS p ≈ 0.01–0.02

--- 

## Current Catalog: 8 Active Systems

| System | Survey | z_l | z_s | Band | Notes |
|--------|--------|-----|-----|------|-------|
| MG0414+0534 | MG | 0.96 | 2.64 | 8.4 GHz | Classic radio quad |
| B0128+437 | CLASS | 0.74 | 3.13 | 5 GHz | |
| B0712+472 | CLASS | 0.41 | 1.34 | 5 GHz | |
| B1422+231 | CLASS | 0.34 | 3.62 | 8.4 GHz | |
| B1608+656 | CLASS | 0.63 | 1.39 | 8.4 GHz | |
| B1933+503 | CLASS | 0.76 | 2.62 | 8.4 GHz | Inner core quad only |
| B2045+265 | CLASS | 0.87 | 1.28 | 8.5 GHz | |
| PG1115+080 | optical | 0.31 | 1.72 | F814W | Optical only |

**Placeholder (inactive):** B1555+375 (redshifts unknown)

---

## Our Simulation Approach

For each draw:
1. Draw source redshift from P(z_s) (CLASS source redshift distribution)
2. Draw lens redshift from P(z_l) (CLASS lens redshift distribution) 
3. Draw lens mass (velocity dispersion) from velocity function
4. Build smooth lens: SIE + external shear, with stellar mass from Chabrier IMF + 2D Sérsic profile
5. Randomly orient the lens
6. Place source in the caustic region to produce 4 images
7. Add either LOS-only or LOS+substructure perturbations
8. Solve lens equation with lenstronomy
9. Compute R_fold and R_cusp from output fluxes
10. Add measurement noise (Gaussian, σ = flux/SNR)
11. Compare simulated distribution to 8 observed systems via KS test

---

## Key Results

### Main Result (5000 draws, ~1900 sim quads, 7 observed — BEFORE catalog update)

| Statistic | KS p-value | Interpretation |
|-----------|-----------|----------------|
| R_fold (without LOS) | 0.0526 | Borderline — smooth model not strongly ruled out |
| R_fold (with LOS) | 0.0145 | Significant — smooth + LOS model inconsistent with data |
| R_cusp (with LOS) | 0.0194 | Significant — independent confirmation |

### After Catalog Update (2000 draws, ~743 sim quads, 8 observed)

| Statistic | KS p-value |
|-----------|-----------|
| R_fold | 0.0496 |
| R_cusp | 0.0029 |

PG1115+080 has R_fold ≈ 0.194 (moderate) but R_cusp ≈ 0.819 (highly anomalous), driving the cusp significance.

### f_sub Posterior
- Flat across f_sub ∈ [0, 0.02]
- LOS alone produces the full anomaly
- Subhalos on top of LOS don't add extra signal
- Implication: f_sub is not constraining with current methodology/data

### Computational Performance
- solve_macro_only: ~14ms/call (with 47 TNFW profiles)
- Full solve: ~1.9s/call (47 TNFW profiles)
- Parallel mode: ~17-19 draws/s on 8 cores
- 5000 draws → ~5 minutes
- --quick mode (50 draws) → ~30s

---

## Code Conventions & Constraints

- **Windows / Python 3.12**: No `fork` — only `spawn` for multiprocessing
- **lenscat**: Must install with `--no-deps` (ligo.skymap broken on Python ≥3.12)
- **lenstronomy profiles**: CONVERGENCE and SHEAR profiles use `ra_0`/`dec_0` kwargs (not `center_x`/`center_y`)
- **Checkpoint system**: `outputs/inference_result.npz` stores incremental results; use `--force` to invalidate on config change
- **CLI flags**: `--quick` (50 draws, parallel, auto-force), `--snr N` (config override), `--n-systems N` (draw override)
- **Git**: All results committed; .planning/ excluded from PR branches

---

## What Should We Do Next?

Please analyze our project status and recommend the next steps. Consider:

1. **Statistical robustness** — Our main result (KS p ≈ 0.01-0.05) hinges on only 8 observed systems. Is this sample size adequate? Should we add more observed systems?
2. **f_sub constraints** — The flat posterior is a problem. Can we modify the analysis (e.g., joint R_fold+R_cusp likelihood, Bayesian hierarchical model) to better constrain f_sub?
3. **Optical vs radio mixing** — PG1115+080 is optical (F814W). Does mixing optical/radio make physical sense? The optical emission region is larger than radio, potentially less sensitive to substructure but more affected by microlensing.
4. **Missing observed systems** — Are there additional CLASS/MG/JVAS radio quads with published fluxes we should add? Known candidates: B0739+366, B2108+213 (but these may be doubles, not quads).
5. **Methodology improvements** — Our smooth lens model (SIE+shear) is relatively simple. How much would more sophisticated models affect the results?
6. **Publication readiness** — What would it take to make this publishable? What plots, tables, and analysis are missing?
7. **Code quality** — Are there specific code improvements, unit tests, or documentation gaps to address before publication?

Please give specific, actionable recommendations ordered by priority.
