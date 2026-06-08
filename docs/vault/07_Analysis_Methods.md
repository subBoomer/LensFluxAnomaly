# Analysis Methods

## Primary: Anderson-Darling Test
The Anderson-Darling test is the preferred comparison statistic as specified by the original project design. It is more sensitive than KS to tail differences, which is where the anomaly signal manifests.

## Secondary: KS Test
Two-sample KS test provides a complementary comparison. Less sensitive to tail differences but more widely understood.

## Tail Ratio
T = P_obs(R_min > 0.2) / P_sim(R_min > 0.2)

The threshold of R_min > 0.2 was specified in the original project design as the definition of the "high tail." This quantifies excess anomaly frequency in the observed sample relative to LCDM.

## Sub-Sample Analysis
Results are computed separately for:
- **Radio-only**: 7 CLASS/MG radio systems (redshifts known, no microlensing)
- **Optical-only**: 10 CASTLES F814W systems (microlensing possible)
- **Combined**: 15 unique systems (deduplicated)

## Comparison with Full Ray-Tracing
The R_min analysis is complemented by a separate analysis using full lenstronomy ray-tracing with R_fold and R_cusp statistics (`run_inference.py`), providing a cross-check of the anomaly signal.

## Code
- `run_rmin_analysis.py` — main entry point with --mode flag
- `run_inference.py` — lenstronomy-based cross-check
- `run_decomposition.py` — diagnostic table for f_sub sensitivity
