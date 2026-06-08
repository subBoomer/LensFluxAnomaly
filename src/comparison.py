import numpy as np
from scipy.stats import ks_2samp, anderson_ksamp
from scipy.stats import wasserstein_distance

def compare(R_obs: np.ndarray, R_sim: np.ndarray, n_bootstrap: int=1000) -> dict:
    ks_stat, ks_p = ks_2samp(R_sim, R_obs, method='exact')
    wass_dist = wasserstein_distance(R_sim, R_obs)
    try:
        ad_result = anderson_ksamp([R_sim, R_obs], variant='midrank')
        ad_stat = ad_result.statistic
        ad_sig = ad_result.critical_values[2]
        ad_p = ad_result.significance_level
    except Exception:
        ad_stat = -1.0
        ad_sig = -1.0
        ad_p = -1.0
    combined = np.concatenate([R_sim, R_obs])
    n_sim = len(R_sim)
    n_obs = len(R_obs)
    rng = np.random.default_rng(42)
    ks_bootstrap = []
    for _ in range(n_bootstrap):
        boot = rng.choice(combined, size=n_sim + n_obs, replace=True)
        ks_bootstrap.append(ks_2samp(boot[:n_sim], boot[n_sim:], method='exact').statistic)
    ks_ci = [float(np.percentile(ks_bootstrap, 2.5)), float(np.percentile(ks_bootstrap, 97.5))]
    return {'n_observed': int(n_obs), 'n_simulated': int(n_sim), 'ks_statistic': float(ks_stat), 'ks_p_value': float(ks_p), 'ks_bootstrap_ci': ks_ci, 'ad_statistic': float(ad_stat), 'ad_significance_level': float(ad_p), 'wasserstein_distance': float(wass_dist)}