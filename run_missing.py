import sys, numpy as np, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from data.radio_quads import build_observed_catalog as get_obs
from src.statistic import compute_rfold
from src.comparison import compare
from run_fsub_grid import prepare_lenses, run_config

outputs_dir = Path(__file__).parent / 'outputs'
samples = prepare_lenses(500, 42)

R_obs = []
for s in get_obs():
    r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
        np.array(s['parity']), np.array(s['fluxes_mjy']))
    if r >= 0:
        R_obs.append(r)
R_obs = np.array(R_obs)
print(f'Observed: {len(R_obs)}')

base = {'substructure': True, 'los': True, 'selection': True, 'noise': True}
for fs in [0.010, 0.020]:
    ckpt = outputs_dir / f'fsub_{fs:.3f}.npy'
    if ckpt.exists():
        print(f'{fs:.3f} cached, skip')
        continue
    cfg_run = {**base, 'f_sub': fs}
    rng = np.random.default_rng(42 + int(fs * 1000))
    t0 = time.perf_counter()
    R_sim = run_config(samples, cfg_run, rng, name=f'f_sub={fs:.3f}')
    np.save(ckpt, R_sim)
    elapsed = time.perf_counter() - t0
    res = compare(R_obs, R_sim, n_bootstrap=1000)
    print(f'f_sub={fs:.3f}: N={len(R_sim)}, KS_p={res["ks_p_value"]:.4f}, Wass={res["wasserstein_distance"]:.4f}, {elapsed:.0f}s')

print('=== All done ===')
