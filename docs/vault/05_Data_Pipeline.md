# Data Pipeline

## Flow
1. **Catalog ingestion** — compile quad systems from radio and optical catalogs
2. **R_min computation** — compute statistic per system
3. **LCDM simulation** — generate null distribution
4. **Statistical comparison** — AD test, KS test, tail ratio

## Implementation
- `src/catalog_utils.py` — unified catalog builder (merges radio + optical)
- `src/compute_rmin.py` — R_min computation from image positions and fluxes
- `src/perturbation_model.py` — simple subhalo perturbation model (no ray tracing)
- `run_rmin_analysis.py` — end-to-end pipeline

## Catalog Sources
- `data/radio_quads.py` — 7 CLASS/MG radio quads
- `data/curated_quads.py` — 10 CASTLES/HST optical quads
- Deduplication: overlapping systems (MG0414, B1422) kept as radio versions

## Sample Selection
- All quads with 4 identified images
- Valid fluxes in a consistent band per system
- Known image positions
- Placeholder systems excluded (B1555 with unknown redshifts)

## R_Min Computation
1. Compute lens centroid as geometric mean of image positions
2. Compute radial distances r_i = sqrt((x_i - x_c)^2 + (y_i - y_c)^2)
3. Normalize fluxes per system: F_i = F_raw_i / sum(F_raw)
4. Find image pairs with |r_i - r_j| < delta_r (default 0.2 arcsec)
5. Compute |F_i - F_j| / (F_i + F_j) for each pair
6. Take the minimum value across all pairs
