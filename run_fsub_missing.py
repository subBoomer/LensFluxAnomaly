import sys, numpy as np, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from run_fsub_grid import prepare_lenses, run_config

def main():
    outputs_dir = Path(__file__).parent / 'outputs'
    from data.radio_quads import build_observed_catalog as get_obs
    from src.statistic import compute_rfold
    systems = get_obs()
    R_obs = []
    for s in systems:
        r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
            np.array(s['parity']), np.array(s['fluxes_mjy']))
        if r >= 0:
            R_obs.append(r)
    R_obs = np.array(R_obs)
    print(f'Observed: {len(R_obs)} systems')

    samples = prepare_lenses(500, 42)
    base_config = {'substructure': True, 'los': True, 'selection': True, 'noise': True}

    for fs in [0.010, 0.020]:
        ckpt = outputs_dir / f'fsub_{fs:.3f}.npy'
        if ckpt.exists():
            print(f'{fs:.3f} already cached, skipping')
            continue
        cfg_run = {**base_config, 'f_sub': fs}
        rng = np.random.default_rng(42 + int(fs * 1000))
        t0 = time.perf_counter()
        R_sim = run_config(samples, cfg_run, rng, name=f'f_sub={fs:.3f}')
        elapsed = time.perf_counter() - t0
        np.save(ckpt, R_sim)
        from src.comparison import compare
        res = compare(R_obs, R_sim, n_bootstrap=1000)
        print(f'  f_sub={fs:.3f}: N={len(R_sim)}, KS_p={res["ks_p_value"]:.4f}, Wass={res["wasserstein_distance"]:.4f}, {elapsed:.0f}s')

    print()
    print('=== Final posterior ===')
    from src.inference import run_inference
    all_ckpts = sorted(outputs_dir.glob('fsub_*.npy'))
    f_grid = []
    R_grid = []
    for ckpt in all_ckpts:
        fs = float(ckpt.stem.replace('fsub_', ''))
        R = np.load(ckpt)
        print(f'  f_sub={fs:.4f}: N={len(R)}, mean={np.mean(R):.4f}')
        if len(R) > 0:
            f_grid.append(fs)
            R_grid.append(R)
    post = run_inference(f_grid, R_grid, R_obs, prior='uniform')
    print(f'\nBest-fit : {post["best_fit"]:.4f}')
    print(f'Median   : {post["median"]:.4f}')
    print(f'68% CI   : [{post["ci_68_low"]:.4f}, {post["ci_68_high"]:.4f}]')
    print(f'Mean±std : {post["mean"]:.4f} ± {post["std"]:.4f}')
    np.savez(outputs_dir / 'fsub_posterior.npz',
        f_sub_vals=post['f_sub_vals'], posterior=post['posterior'],
        best_fit=np.array([post['best_fit']]), median=np.array([post['median']]),
        ci_68_low=np.array([post['ci_68_low']]), ci_68_high=np.array([post['ci_68_high']])
    )

if __name__ == '__main__':
    main()
