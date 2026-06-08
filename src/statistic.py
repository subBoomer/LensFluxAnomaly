import numpy as np


def compute_rfold(theta_x: np.ndarray, theta_y: np.ndarray,
                  mu: np.ndarray, F_obs: np.ndarray) -> float:
    """Classical fold relation: R_fold = |F_a - F_b| / (F_a + F_b)

    The fold pair is the two minimum (positive-parity) images.
    """
    n = len(theta_x)
    if n < 4:
        return -1.0
    parity = np.array([1 if m > 0 else -1 for m in mu])
    pos = np.where(parity == 1)[0]
    if len(pos) < 2:
        return -1.0
    i, j = int(pos[0]), int(pos[1])
    denom = F_obs[i] + F_obs[j]
    if denom <= 0:
        return -1.0
    return float(abs(F_obs[i] - F_obs[j]) / denom)


def compute_cusp_statistic(theta_x: np.ndarray, theta_y: np.ndarray,
                           mu: np.ndarray,
                           F_obs: np.ndarray | None = None) -> float:
    """Cusp statistic: R_cusp = |sum(signed)| / sum(|signed|)

    For simulated data (mu = signed magnifications, F_obs=None):
        R_cusp = |mu_a + mu_b + mu_c| / (|mu_a| + |mu_b| + |mu_c|)

    For observed data (mu = parity array, F_obs = fluxes):
        R_cusp = |p_a*F_a + p_b*F_b + p_c*F_c| / (F_a + F_b + F_c)

    Cusp triplet: the three closest images (standard literature criterion).

    NOTE: The cusp result is sensitive to triplet identification.
    See PLAN.md for stability analysis showing the rejection direction
    flips depending on whether close same-parity pairs are included.
    """
    n = len(theta_x)
    if n < 4:
        return -1.0
    min_sum = float('inf')
    triplet = None
    for i in range(n):
        for j in range(i + 1, n):
            d_ij = np.hypot(theta_x[i] - theta_x[j], theta_y[i] - theta_y[j])
            for k in range(j + 1, n):
                d_sum = d_ij + np.hypot(theta_x[i] - theta_x[k],
                                        theta_y[i] - theta_y[k]) + \
                        np.hypot(theta_x[j] - theta_x[k],
                                 theta_y[j] - theta_y[k])
                if d_sum < min_sum:
                    min_sum = d_sum
                    triplet = (i, j, k)
    i, j, k = triplet
    if F_obs is not None:
        vals = np.array([mu[i] * F_obs[i], mu[j] * F_obs[j],
                         mu[k] * F_obs[k]])
    else:
        vals = np.array([mu[i], mu[j], mu[k]])
    denom = np.sum(np.abs(vals))
    if denom <= 0:
        return -1.0
    return float(abs(np.sum(vals)) / denom)
