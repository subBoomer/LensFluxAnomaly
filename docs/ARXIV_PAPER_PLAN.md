# arXiv Paper Plan

## What the Paper Says

A short (~8 page) paper reporting that the R_fold flux-ratio distribution in 8 CLASS radio quads is inconsistent with CDM+substructure predictions at p < 0.001, and that this result is robust to selection bias, model choice, and jackknife, but is not yet a cosmological constraint due to sample size.

## Authorship

Your Name (sole author, or "LensFluxAnomaly Collaboration" if preferred)

## Structure

### 1. Introduction (1.5 pages)
- Flux-ratio anomalies as substructure probes (Dalal & Kochanek 2002)
- Previous work: CLASS sample, CASTLES, Nierenberg+2017, Hsueh+2020
- Gap: most studies fit individual lenses, we build a forward model for the population
- Our approach: three statistics, two models, one catalog, honest error budget

### 2. Data (1.5 pages)
- CLASS radio sample (N=8): table with positions, fluxes, parity, references
- Selection function: CLASS selected on total flux + spectral index, NOT flux ratios
- Selection bias test: 49,980 CDM realizations → 0% bias (confirmed)
- Optical CASTLES sample (N=13, diagnostic only)

### 3. Methods (2 pages)
- Three statistics: R_min, R_fold (primary), R_cusp
- Simple perturbation model: dN/dm ~ m^{-1.9}, epsilon kernel, 50k realizations
- Lenstronomy pipeline: SIE+TNFW+LOS, 743 realizations, f_sub=0.005
- WDM: 3/5/7 keV, mass function suppression
- Statistical tests: AD (primary), KS (secondary), PPC (calibrated), jackknife

### 4. Results (1.5 pages)
- R_fold: observed mean 0.436 vs CDM 0.165 (2.6x), AD p=0.001, PPC p<10^{-4}
- Lenstronomy strengthens all signals (tighter null)
- WDM makes anomaly worse
- Jackknife: signal survives removal of any single system
- Per-system table with R_min, R_fold, R_cusp

### 5. Discussion (1 page)
- What we can and cannot conclude
- N=8 is insufficient for f_sub constraints
- Future: VLASS, Euclid, Roman will expand samples
- R_fold recommended as primary statistic for future work

### 6. Conclusion (0.5 pages)
- Observed R_fold distribution is inconsistent with minimal CDM at high significance
- Selection, model, and outlier-robustness checks confirm the result
- Larger samples needed for cosmological interpretation

## Figures Needed
1. **R_fold distribution**: histogram, observed vs CDM vs lenstronomy (primary figure)
2. **R_min distribution**: same for combined sample (secondary)
3. **PPC distribution**: observed mean vs 10,000 CDM samples (calibration proof)
4. **Jackknife**: AD statistic when each system is removed
5. **WDM comparison**: R_min means for CDM vs 3/5/7 keV WDM

All figures already generated in `outputs/`.

## Timeline

| Step | Time |
|------|------|
| Write introduction | 1 day |
| Write data section | 1 day |
| Write methods | 2 days |
| Write results (with figures) | 1 day |
| Write discussion + conclusion | 1 day |
| Compile references, format, abstract | 1 day |
| Review + finalize | 1-2 days |
| **Total** | **~8-10 days** |

## What I'd need from you
- Your name + affiliation for the author line
- Which figure style you prefer (color scheme, etc.)
- Any specific references you want cited
- Final review of the draft

## Submission
- Venue: arXiv (astro-ph.GA)
- No journal submission needed
- Can update later if you want to submit to a journal
