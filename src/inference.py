import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde


def build_emulator(f_sub_grid, R_grid, bw_method='scott'):
    grid = np.asarray(f_sub_grid, dtype=float)
    n = len(grid)
    log_pdfs = []
    for R_sim in R_grid:
        R_sim = np.asarray(R_sim, dtype=float)
        if len(R_sim) < 5:
            log_pdfs.append(None)
            continue
        kde = gaussian_kde(R_sim, bw_method=bw_method)
        log_pdfs.append(kde)

    def log_likelihood(f_sub, R_obs):
        idx = np.searchsorted(grid, f_sub)
        if idx <= 0:
            kde_left = log_pdfs[0]
            if kde_left is None:
                return -np.inf
            return float(np.sum(np.log(np.maximum(kde_left.evaluate(R_obs), 1e-300))))
        if idx >= n:
            kde_right = log_pdfs[-1]
            if kde_right is None:
                return -np.inf
            return float(np.sum(np.log(np.maximum(kde_right.evaluate(R_obs), 1e-300))))
        left_kde = log_pdfs[idx - 1]
        right_kde = log_pdfs[idx]
        if left_kde is None or right_kde is None:
            return -np.inf
        f_left = grid[idx - 1]
        f_right = grid[idx]
        w = (f_sub - f_left) / (f_right - f_left) if f_right > f_left else 0.5
        log_left = np.sum(np.log(np.maximum(left_kde.evaluate(R_obs), 1e-300)))
        log_right = np.sum(np.log(np.maximum(right_kde.evaluate(R_obs), 1e-300)))
        return float((1.0 - w) * log_left + w * log_right)

    return log_likelihood, grid


def posterior_grid(f_sub_grid, R_obs, log_likelihood, prior='uniform', n_grid=1000):
    f_min = float(np.min(f_sub_grid))
    f_max = float(np.max(f_sub_grid))
    f_sub_vals = np.linspace(f_min, f_max, n_grid)
    log_post = np.zeros(n_grid)
    for i, fs in enumerate(f_sub_vals):
        log_l = log_likelihood(fs, R_obs)
        if prior == 'log_uniform':
            log_prior = -np.log(fs) if fs > 0 else -np.inf
        else:
            log_prior = 0.0
        log_post[i] = log_l + log_prior
    log_post -= np.max(log_post)
    post = np.exp(log_post)
    post /= np.trapezoid(post, f_sub_vals)
    mean = np.trapezoid(post * f_sub_vals, f_sub_vals)
    var = np.trapezoid(post * (f_sub_vals - mean) ** 2, f_sub_vals)
    cdf = np.cumsum(post) * (f_sub_vals[1] - f_sub_vals[0])
    idx_median = np.searchsorted(cdf, 0.5)
    median = f_sub_vals[idx_median]
    idx_low = np.searchsorted(cdf, 0.16)
    idx_high = np.searchsorted(cdf, 0.84)
    ci_lo = f_sub_vals[idx_low]
    ci_hi = f_sub_vals[idx_high]
    idx_best = np.argmax(post)
    best = f_sub_vals[idx_best]
    return {
        'f_sub_vals': f_sub_vals,
        'posterior': post,
        'mean': float(mean),
        'median': float(median),
        'best_fit': float(best),
        'ci_68_low': float(ci_lo),
        'ci_68_high': float(ci_hi),
        'std': float(np.sqrt(var)),
    }


def run_inference(f_sub_grid, R_grid, R_obs, prior='uniform', bw_method='scott'):
    log_likelihood, grid = build_emulator(f_sub_grid, R_grid, bw_method=bw_method)
    result = posterior_grid(grid, R_obs, log_likelihood, prior=prior)
    result['f_sub_grid'] = list(f_sub_grid)
    result['log_likelihood'] = log_likelihood
    return result
