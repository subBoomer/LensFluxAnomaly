---
tags:
  - analysis
  - methods
  - lensfluxanomaly
---
# Analysis Methods

See [[06_Statistical_Methods]] for the statistical testing framework.

## Entry-Point Scripts
| Script | Function | Approx Time |
|--------|----------|-------------|
| `run_rmin_analysis.py` | R_min for all/radio/optical | ~30 sec |
| `run_rfold_rcusp_analysis.py` | R_fold/R_cusp for radio | ~30 sec |
| `run_phase_f.py` | WDM comparison | ~3 min |
| `run_validate_model.py --quick` | lenstronomy cross-validation | ~30 sec |
| `run_validate_model.py` | Full lenstronomy (2000) | ~5-10 min |
| `run_jackknife.py` | Jackknife R_min | ~2 min |
| `run_jackknife_rfrc.py` | Jackknife R_fold/R_cusp | ~2 min |
| `run_outlier_analysis.py` | Individual system diagnostics | ~30 sec |

All scripts use `pathlib.Path(__file__).parent` for path resolution and can be run from any working directory.

## Verification
The four critical bugs (global rng, file mutation, hash seeding, float rounding) have been fixed as of the final analysis run.
