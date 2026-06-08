import numpy as np
from pathlib import Path
from data.radio_quads import build_observed_catalog
from src.statistic import compute_rfold
from src.comparison import compare

ckpt = np.load('outputs/checkpoint.npz')
R_sim = list(ckpt['R_sim'])

observed_systems = build_observed_catalog()
R_obs = []
for sys in observed_systems:
    parity = np.array(sys['parity'])
    r = compute_rfold(
        np.array(sys['x_arcsec']), np.array(sys['y_arcsec']),
        parity, np.array(sys['fluxes_mjy']),
    )
    if r >= 0:
        R_obs.append(r)
R_obs = np.array(R_obs)
R_sim = np.array(R_sim)
result = compare(R_obs, R_sim, 1000)

save_dict = {
    'R_obs': R_obs, 'R_sim': R_sim,
    'n_observed': np.array([result['n_observed']]),
    'n_simulated': np.array([result['n_simulated']]),
    'ks_statistic': np.array([result['ks_statistic']]),
    'ks_p_value': np.array([result['ks_p_value']]),
    'wasserstein_distance': np.array([result['wasserstein_distance']]),
}
np.savez('outputs/inference_result.npz', **save_dict)
try:
    Path('outputs/checkpoint.npz').unlink()
except PermissionError:
    pass

print(f'Observed: {len(R_obs)} systems')
for s, r in zip(observed_systems, R_obs):
    print(f'  {s["name"]}: R={r:.4f}  parity={s["parity"]}')
print(f'Simulated: {len(R_sim)}')
print(f'KS: stat={result["ks_statistic"]:.4f}, p={result["ks_p_value"]:.4f}')
print(f'Wasserstein: {result["wasserstein_distance"]:.4f}')
