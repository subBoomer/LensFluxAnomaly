# Validation + Conclusion Plan

## Objective
Answer: **"Is there a ΛCDM flux-ratio anomaly, quantified, with caveats?"**

## Phases

### Phase A: Expand Observed Catalog
Add ~15-20 more optical quads from SQLS (Inada+2012, J/AJ/143/119) and HSC (Jaelani+2020) to bring the sample from N=15 to N=30+.

**Files**: `data/curated_quads.py` (add entries in same format)

### Phase B: Cross-Validate Models
Compute R_min from the full lenstronomy pipeline (SIE+TNFW+LOS+SNR) — our most realistic ΛCDM model — and compare:
1. Lenstronomy R_min vs observed R_min
2. Lenstronomy R_min vs simple perturbation model R_min

**File**: `run_validate_model.py`

### Phase C: Outlier + Cross-Statistic Table
For all systems:
- R_min value + tail flag
- R_fold value (radio systems only)
- Radio vs optical R_min comparison for overlapping systems (MG0414, B1422)

**File**: `run_outlier_analysis.py`

### Phase D: Jackknife Stability
Remove each system one at a time, re-run AD test + tail ratio. Track stability.

**File**: `run_outlier_analysis.py` (same script)

### Phase E: Conclusion Writeup
One-page answer: detection significance, driving systems, caveats, path to definitive answer.

**File**: `docs/CONCLUSION.md`
