# LensFluxAnomaly — Project Plan

## 1. Current Status

| Metric | Value |
|--------|-------|
| Simulated draws | 2000 |
| Valid quads (post-selection) | 806 |
| Observed quads | 7 (CLASS/MG) |
| KS p-value | **0.056** (borderline) |
| Wasserstein distance | 0.205 |
| Ablation configs | 6/6 complete |
| f_sub response grid | 4 points (0, 0.002, 0.005, 0.01) |

**Interpretation**: KS p=0.056 is tantalizing but not significant at the 5% threshold. The result is consistent with substructure at `f_sub ≈ 0.005` but does not rule out smooth models with >5% probability.

---

## 2. Performance Bottlenecks

### Current Timing (2000 draws → ~25 min run time)

| Component | Per-call | Per-2000-draws | % of total |
|-----------|----------|----------------|------------|
| `solve_macro_only` (30 attempts) | ~55ms ea → ~1.0-2.5s/draw | ~2000-5000s | ~60% |
| Full `solve` (with subhalos) | ~1.9s | ~1500s (for ~800 valid) | ~30% |
| Population, build, selection, noise | negligible | ~10s | <1% |

### Primary bottleneck: `solve_macro_only` creates a fresh `LensModel` on every call
File: `src/lens_model.py:87-93` — each call instantiates:
```python
fast_solver = LensEquationSolver(LensModel(['SIE', 'SHEAR']))
```
This is unnecessary — the model doesn't change during caustic sampling.

### Secondary bottleneck: Full solve with 49 TNFW subhalos
`LensModel` with 51 profiles (SIE + SHEAR + 49x TNFW) is slow to ray-trace through.

### No parallelization
The 2000-draw loop in `run_inference.py:85` is strictly sequential.

### Quick Wins (implement first, measure improvement)

1. **Cache the macro solver** in `solve_macro_only` — store `self._macro_solver` and `self._macro_model` during `build()`, use them instead of creating fresh ones:
   ```
   Current: 55ms/call × 30 att × 2000 draws ≈ 3300s
   After fix: same (caching doesn't speed the solver, only removes instantiation overhead ~2-3ms)
   ```
   Actually caching gives minimal gain. The real cost is in `find_bright_image`.

2. **Reduce `num_random` and `search_window`** in `solve_macro_only`:
   - Current: `num_random=5, search_window=5` → ~79ms
   - Lower: `num_random=3, search_window=3` → ~57ms (~28% faster)
   - Lower: `num_random=2, search_window=2` → ~51ms (~35% faster)
   
3. **Reduce `max_attempts`** in `sample_source_in_caustic`:
   - Current: 30 attempts — tested each via `solve_macro_only`
   - Most quads found in 1-10 attempts. A truncated search (max 15-20) may lose <5% of quads.

### Big Wins (larger refactors)

4. **Multiprocessing** — 2000 draws are completely independent:
   - Use `concurrent.futures.ProcessPoolExecutor` on N workers
   - ~4x speedup on 8-core desktop (less if solver is GPU-accelerated)
   - Challenge: pickling `LensModel` objects; use per-worker initialization

5. **Shared caustic-sampled lens library** — `run_ablation.py` already does this: pre-sample all quads once, then evaluate different model configs on the same set. Port this approach to the main pipeline.

6. **Subhalo count control** — Currently `SubhaloPopulation.realise()` samples from Poisson(λ). For most systems λ≈50, yielding ~49 subhalos on average. Each TNFW profile adds solver cost. Investigate if fewer profiles with scaled masses preserves the R_fold distribution.

---

## 3. Scientific Refinements

### 3.1 Line-of-Sight (LOS) — HIGH priority
**Where**: `src/substructure.py:69-77`
**What**: `LOSPopulation` is a stub returning `kappa=0, gamma_1=0, gamma_2=0`.
**Why import**: LOS halos contribute ~30% of the total substructure convergence in the literature (Despali+2022, Lazar+2021). Ignoring them systematically biases the comparison.
**How**: Implement using `colossus` mass function integrated along the line of sight, project halos as convergence + shear pertubations.

### 3.2 Constrain f_sub via MCMC — HIGH priority
**Where**: `src/comparison.py` → new file `src/inference.py`
**What**: Current grid search (4 f_sub values) is coarse. Use Approximate Bayesian Computation (ABC) or emulation to likelihood-free infer f_sub with uncertainties.
**Why**: KS p=0.056 is inconclusive. A proper posterior on f_sub would show whether it's truly >0.
**How**: Run simulations at [0.0, 0.001, 0.002, ..., 0.020], interpolate R_fold distribution, compute likelihood for observed data, sample posterior.

### 3.3 Increase simulated sample size — MEDIUM priority
**What**: Current 806 valid quads from 2000 draws. Target 4000+ valid quads (config `n_systems: 10000`).
**Why**: Larger simulation sample reduces Monte Carlo noise in the KS test, tightening the constraint.
**How**: Increase `n_systems` in `config.yaml` after implementing performance improvements.

### 3.4 Cusp statistic in main pipeline — LOW priority
**Where**: `src/statistic.py:24-66` (exists but unused in main pipeline)
**Why**: Cusp statistic provides an independent test of the substructure hypothesis. If both R_fold and R_cusp show the same trend, confidence increases.

### 3.5 Add more observed systems — MEDIUM priority
**Where**: `data/radio_quads.py` (current catalog of 7 active + 2 placeholders)
**Considerations**:
- **B1933+503** (placeholder): Has astrometry but positions all zeros; needs proper image coordinates
- **PG1115+080** (placeholder): Optical-only, not radio; higher systematic risk
- **New candidates**: Search CLASS survey archives for additional quads
- **Optical quads**: `data/curated_quads.py` has 10 optical systems from CASTLES — could be used with caveats (different selection function, band)

### 3.6 Systematic error analysis — MEDIUM priority
- **Mass slope sensitivity**: Test `mass_slope ∈ [-1.7, -1.9, -2.1]` (`config.yaml` line 33)
- **Concentration model**: Test `duffy08`, `ishiyama21`, `diemer15` (`config.yaml` line 35)
- **SNR sensitivity**: Test `radio_snr ∈ [20, 50, 100]`
- **Parity ambiguity**: Some systems have disputed parity assignments (see `data/radio_quads.py` parity_source comments)

### 3.7 Better selection function — LOW priority
Current `passes_selection` uses a simple flux cut + quad detection probability. CLASS survey selection is more complex (resolution, spectral index, etc.). Could port to a more realistic survey simulator.

---

## 4. Code Quality & Engineering

### 4.1 `--force` flag for re-runs
**Where**: `run_inference.py` (lines 68-70 reference it but it's not implemented)
**Why**: Cannot re-run without manually deleting `inference_result.npz`

### 4.2 Config validation
**Where**: `config.yaml` → new `src/config.py`
**Why**: Silent fallbacks for missing keys (e.g., `cfg.get('checkpoint_interval', 100)`) can hide bugs

### 4.3 Logging
**Current**: Bare `print()` statements
**Target**: Use `logging` module with timestamps, log levels, output to file when needed

### 4.4 Tests
**Where**: New `tests/` directory
**Why**: Zero tests currently. Critical for regression when optimizing.
- Unit tests for `compute_rfold` (known test cases)
- Unit tests for `LensPopulation.sample()` (distribution shapes)
- Integration test for one full pipeline iteration

### 4.5 Ablation visualization
**Where**: `src/plot_results.py` or new `src/plot_ablation.py`
**What**: Compare R_fold distributions across all 6 ablation configs on one figure
**Why**: Current output only has `lens_anomaly_plot.png` — no comparison of model complexity

---

## 5. Implementation Roadmap

### Phase A — Quick Fixes (1-2 days)
- [ ] Implement `--force` flag in `run_inference.py`
- [ ] Cache macro solver in `solve_macro_only` (avoid re-instantiation)
- [ ] Add config validation
- [ ] Run with `n_systems=5000` to get ~2000 valid quads

### Phase B — Performance (2-3 days)
- [ ] Multiprocess the draw loop in `run_inference.py`
- [ ] Benchmark speedup and tune chunk size
- [ ] Test subhalo count sensitivity (fewer but scaled masses)

### Phase C — Scientific Depth (1-2 weeks)
- [ ] Implement LOS population (real halos, not zeros)
- [ ] f_sub MCMC inference
- [ ] Add cusp statistic to main pipeline
- [ ] Create ablation comparison plots
- [ ] Systematic tests (mass slope, concentration, SNR)

### Phase D — Robustness (ongoing)
- [ ] Add observed systems (fix placeholders, search for more)
- [ ] Write unit tests
- [ ] Logging improvements
- [ ] Comment parity provenance for all observed systems

---

## 6. Data Flow

```
config.yaml
    ↓
LensPopulation.sample() ───────────────────── 2000× draws
    ↓                                          (independent, parallelizable)
MacroLens.build() → theta_E
    ↓
sample_source_in_caustic() ←────────────────  up to 30 attempts via solve_macro_only
    ↓                                          (major bottleneck)
passes_selection() ────────────────────────── flux cut, resolution, quad detection
    ↓
SubhaloPopulation.realise() ───────────────── Poisson(λ=50) TNFW profiles
    ↓
MacroLens.add_substructure()
    ↓
MacroLens.solve() ────────────────────────── full lens equation with 51 profiles
    ↓                                          (second bottleneck)
RadioNoise.apply()
    ↓
compute_rfold()
    ↓
[R_sim] ──→ compare(R_obs, R_sim) ────────── KS + Wasserstein + bootstrap
```

---

## 7. File Reference

| File | Purpose | Needs work? |
|------|---------|-------------|
| `run_inference.py` | Main pipeline entry point | --force flag, multiprocessing, cusp stat |
| `run_ablation.py` | Ablation matrix + f_sub response | Already cached results, may need re-run |
| `finalize.py` | Quick finalize from checkpoint | Redundant with run_inference (remove?) |
| `src/lens_model.py` | SIE+SHEAR lens + solver | Cache solver, reduce search params |
| `src/substructure.py` | Subhalo + LOS population | **LOS is a stub** — highest scientific gap |
| `src/statistic.py` | R_fold + R_cusp computation | Add cusp to main pipeline |
| `src/comparison.py` | KS + Wasserstein + bootstrap | Add more comparison stats |
| `src/population.py` | Lens population sampling | Distribution tuning |
| `src/selection.py` | CLASS survey selection | Better survey model |
| `src/noise_model.py` | Radio noise (Gaussian) | Non-Gaussian noise? |
| `src/plot_results.py` | Results visualization | Add ablation comparison |
| `data/radio_quads.py` | Observed radio quad catalog | Fix B1933 placeholder, add more |
| `data/curated_quads.py` | Optical quads (CASTLES) | Potential additional data |
| `src/ingest.py` | Data ingestion pipeline | Not used in main pipeline |
| `src/lens_filters.py` | Catalog filtering utilities | Not used in main pipeline |
| `config.yaml` | All configurable parameters | Needs validation |

---

## 8. Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| Initial | Use SIE+SHEAR lens model | Standard in literature for CLASS quads |
| Initial | f_sub = 0.005 for DMO | Springel 2008, Fiacconi+2016 |
| Initial | TNFW for subhalos | Truncated profile prevents infinite mass |
| Initial | 2000 draws × 30 attempts | Balance between compute time and sample size |
| Initial | 7 observed quads | All CLASS/MG radio quads with published parity |
| This plan | Flag LOS as highest gap | Known ~30% effect ignored |
| This plan | Multiprocessing as priority | Quickest path to larger N |
