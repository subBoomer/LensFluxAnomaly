import numpy as np

BETA_MAGNIFICATION = 2.3
FLUX_LIMIT_MJY = 1.0
BEAM_ARCSEC = 0.25
QUAD_DETECT_PROB = 0.7


def passes_selection(theta_x, theta_y, mu, F_source, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    if len(theta_x) < 4:
        return False
    mu_tot = np.sum(np.abs(mu))
    lensed_flux = mu_tot * F_source
    if lensed_flux < FLUX_LIMIT_MJY:
        return False
    min_sep = float('inf')
    n = len(theta_x)
    for i in range(n):
        for j in range(i + 1, n):
            sep = np.hypot(theta_x[i] - theta_x[j], theta_y[i] - theta_y[j])
            min_sep = min(min_sep, sep)
    if n == 4 and min_sep < BEAM_ARCSEC:
        p_resolution = 1.0 / (1.0 + np.exp(-(min_sep - BEAM_ARCSEC) / 0.05))
        if rng.uniform() > p_resolution:
            return False
    return rng.uniform() < QUAD_DETECT_PROB


def selection_weight(theta_x, theta_y, mu, F_source):
    if len(theta_x) < 4:
        return 0.0
    mu_tot = np.sum(np.abs(mu))
    lensed_flux = mu_tot * F_source
    if lensed_flux < FLUX_LIMIT_MJY:
        return 0.0
    n = len(theta_x)
    min_sep = float('inf')
    for i in range(n):
        for j in range(i + 1, n):
            sep = np.hypot(theta_x[i] - theta_x[j], theta_y[i] - theta_y[j])
            min_sep = min(min_sep, sep)
    p_resolution = 1.0 / (1.0 + np.exp(-(min_sep - BEAM_ARCSEC) / 0.05))
    p_quad = QUAD_DETECT_PROB
    return lensed_flux ** BETA_MAGNIFICATION * p_resolution * p_quad
