# Execution Plan

## Step 1: Fix 4 Critical Bugs (~2 hrs)

### 1a. `passes_selection()` — add rng parameter
- `src/selection.py`: add `rng=None` parameter, use `rng.uniform()` instead of `np.random.uniform()`
- `run_inference.py`, `run_validate_model.py`, `run_decomposition.py`, `src/run_ablation.py`: pass local rng

### 1b. `run_systematic.py` — stop mutating source files
- `src/substructure.py`: add `concentration_model` parameter to SubhaloPopulation
- `run_systematic.py`: remove file read/write, pass parameter directly

### 1c. `hash(name)` → `zlib.crc32(name.encode())`
- `src/run_ablation.py`: replace `hash(name)` with deterministic crc32

### 1d. Float rounding in seed derivation
- `run_fsub_grid.py`: `int(round(fs * 1000))` instead of `int(fs * 1000)`

## Step 2: Lenstronomy R_fold/R_cusp (~20 min code + overnight run)

### 2a. Modify `run_validate_model.py`
- Compute R_fold/R_cusp from magnifications in lenstronomy loop
- Save to npz output

### 2b. Run
```bash
python run_validate_model.py --n-sim 2000 --parallel
```

### 2c. Analysis
Compare observed vs lenstronomy R_fold/R_cusp

## Step 3: Build Obsidian Vault (~4-5 hrs)
- Create/update 10 markdown files in docs/vault/

## Step 4: Selection Function Check (~3-4 days)
### 4a. Build CLASS selection model from literature
### 4b. Integrate into forward model
### 4c. Bias measurement simulation
### 4d. Report results

## Step 5: Posterior Predictive Check (~1 day)
### 5a. Implement PPC function
### 5b. Compare with KS/AD p-values
### 5c. Create `run_ppc.py`

## Step 6: README Overhaul (~2 hrs)

## Execution Order
```
Step 1 (all 4 sub-steps in parallel)
    |
    ▼
Step 2a → Step 2b (overnight) → Step 2c
    |
    ▼
Step 5 (PPC, uses outputs from Steps 1-2)
    |
    ├── Step 3 (vault, parallel with Step 4)
    ├── Step 4 (selection, parallel with Step 3)
    └── Step 6 (README, last)
```
