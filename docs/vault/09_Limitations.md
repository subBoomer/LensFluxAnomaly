# Limitations

## Sample Size
With only 15 systems (7 radio + 10 optical with 2 overlap), the statistical power is limited. The tail ratio T = 8.0 is suggestive but the error bars are large. A sample of 50+ quads would provide a definitive test.

## Simple Perturbation Model
The phenomenological kernel (epsilon * m / (r^2 + r_c^2)) is a crude approximation:
- No ray tracing — cannot capture magnification bias or caustic effects
- Equal macro fluxes assumed (F_i_0 = 0.25) — real lenses have intrinsic asymmetries
- No explicit LOS component — line-of-sight structure is absorbed into the kernel
- Calibrated to reproduce small perturbations, not validated against full simulations

## Optical Contamination
Optical systems are subject to:
- Microlensing by stars in the lens galaxy (can cause 10-30% variability)
- Dust extinction along different image paths
- Photometric errors
These effects inflate R_min values, potentially mimicking a substructure signal.

## Radio Sample Heterogeneity
The 7 radio systems come from different surveys (MG, CLASS) at different frequencies (5-8.5 GHz). While radio is less affected by microlensing, variability in AGN cores over observation epochs could contribute.

## Lens Centroid Approximation
The geometric mean of image positions is used as the lens centroid. This can be offset from the true lens center, especially for asymmetric configurations, introducing systematic error in r_i.

## f_Sub Insensitivity
The R_fold/R_cusp statistics (and by extension R_min) are fundamentally limited in their ability to constrain substructure abundance. The macro-model variance dominates the signal, making it impossible to distinguish f_sub values with this approach alone.

## No Blind Analysis
The sample was compiled with knowledge of the science case. A blind analysis would be preferable to avoid confirmation bias.
