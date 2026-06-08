import numpy as np


def compute_rmin(x_arcsec, y_arcsec, fluxes, delta_r=0.2, centroid=None):
    """R_min = min |F_i - F_j| / (F_i + F_j) for pairs with |r_i - r_j| < delta_r.

    Parameters
    ----------
    x_arcsec : array-like
        Image x positions in arcsec.
    y_arcsec : array-like
        Image y positions in arcsec.
    fluxes : array-like
        Image fluxes (any units, will be normalised per system).
    delta_r : float
        Radial pairing threshold in arcsec.
    centroid : tuple or None
        (x_c, y_c) lens centroid. If None, use geometric mean of image positions.

    Returns
    -------
    rmin : float or None
        R_min value, or None if no valid pair found.
    """
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    F = np.asarray(fluxes, dtype=float)
    n = len(x)
    if n < 4:
        return None
    F_norm = F / np.sum(F)
    if centroid is not None:
        x_c, y_c = centroid
    else:
        x_c, y_c = np.mean(x), np.mean(y)
    r = np.hypot(x - x_c, y - y_c)
    best = None
    for i in range(n):
        for j in range(i + 1, n):
            if abs(r[i] - r[j]) < delta_r:
                val = abs(F_norm[i] - F_norm[j]) / (F_norm[i] + F_norm[j])
                if best is None or val < best:
                    best = val
    return best


def compute_rmin_from_catalog(systems, delta_r=0.2):
    """Compute R_min for each system in a catalog.

    Parameters
    ----------
    systems : list of dict
        Each dict must have keys 'x_arcsec', 'y_arcsec', 'fluxes_mjy', 'name'.
    delta_r : float
        Radial pairing threshold in arcsec.

    Returns
    -------
    results : list of (name, rmin) tuples
    """
    results = []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes_mjy'], delta_r)
        results.append((s['name'], rmin))
    return results
