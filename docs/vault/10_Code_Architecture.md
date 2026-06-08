# Code Architecture

## Top-Level Scripts
| Script | Purpose |
|--------|---------|
| `run_rmin_analysis.py` | R_min analysis with simple perturbation model |
| `run_inference.py` | R_fold/R_cusp analysis with full lenstronomy |
| `run_decomposition.py` | Diagnostic: f_sub sensitivity test |
| `run_systematic.py` | Systematic tests (SNR sweep, concentration sweep) |
| `run_fsub_grid.py` | f_sub grid search |

## Source Modules
| Module | Purpose |
|--------|---------|
| `src/compute_rmin.py` | R_min statistic computation |
| `src/perturbation_model.py` | Simple subhalo perturbation model |
| `src/catalog_utils.py` | Unified catalog builder |
| `src/statistic.py` | R_fold and R_cusp statistics |
| `src/lens_model.py` | Macro lens model (lenstronomy SIE) |
| `src/substructure.py` | Subhalo and LOS population |
| `src/population.py` | Lens population sampling |
| `src/selection.py` | Selection function |
| `src/noise_model.py` | Radio noise model |
| `src/comparison.py` | Statistical comparison utilities |
| `src/inference.py` | f_sub inference (KDE + grid) |
| `src/config.py` | Configuration validation |
| `src/ingest.py` | Catalog ingestion (online) |
| `src/lens_filters.py` | Catalog filtering |
| `src/literature_catalog.py` | Literature catalog builder |

## Data Files
| File | Contents |
|------|----------|
| `data/radio_quads.py` | 7 curated radio quad systems |
| `data/curated_quads.py` | 10 CASTLS optical quads |

## Outputs
| File | Contents |
|------|----------|
| `outputs/rmin_comparison.png` | R_min distribution plot |
| `outputs/rmin_radio.png` | Radio-only R_min plot |
| `outputs/rmin_optical.png` | Optical-only R_min plot |
| `outputs/rmin_all.png` | Combined R_min plot |
| `outputs/inference_result.npz` | R_fold/R_cusp inference results |

## Dependencies
- Python 3.12+
- numpy, scipy, matplotlib
- lenstronomy (for ray-tracing pipeline only)
- astropy (for catalog utilities)
