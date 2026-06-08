---
tags:
  - pipeline
  - data
  - lensfluxanomaly
---
# Data Pipeline

## Ingestion
- Radio systems from `data/radio_quads.py` (manually curated from literature)
- Optical systems from `data/curated_quads.py` (CASTLES HST F814W)
- Unified via `src/catalog_utils.py:build_unified_catalog()`
- Deduplication by name (radio takes priority over optical)

## Position Normalization
Image positions and lens centroid positions are taken directly from literature.
No re-centering or transformation applied.

## Flux Normalization
Optical magnitudes converted to linear flux: f = 10^{-0.4 * mag}
All fluxes used in raw units (not normalized to sum=1) for R_min/R_fold/R_cusp computation.

## Statistics Computation
- R_min: `src/compute_rmin.py` — pairs images within δr = 0.2 arcsec
- R_fold: `src/statistic.py` — identifies positive-parity images
- R_cusp: `src/statistic.py` — identifies closest triple via minimum distance sum

## Output
Results saved to `outputs/*.npz` files with standard key names.
Plots saved to `outputs/*.png`.

See [[10_Code_Architecture]] for detailed file listing.
