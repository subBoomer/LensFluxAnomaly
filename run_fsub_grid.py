import sys, numpy as np, yaml, time, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from src.config import validate as validate_config
from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation, LOSPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from src.statistic import compute_rfold
from src.comparison import compare


def prepare_lenses(n_total, seed):
    rng = np.random.default_rng(seed)
    pop = LensPopulation(rng)
    samples = []
    for i in range(n_total):
        theta = pop.sample(rng)
        macro = MacroLens()
        theta_E = macro.build(theta['z_l'], theta['z_s'],
            theta['sigma_v'], theta['q'], theta['phi_lens'],
            theta['gamma_ext'], theta['theta_gamma'])
        cs = theta_E * (1.0 - theta['q']) / (1.0 + theta['q'])
        scale = max(cs * 0.8, 0.01)
        for _ in range(30):
            bx = rng.uniform(-scale, scale)
            by = rng.uniform(-scale, scale)
            result = macro.solve_macro_only(bx, by)
            if result is not None and result['n_images'] == 4:
                samples.append((theta, theta_E, bx, by, result))
                break
        if (i + 1) % 100 == 0:
            print(f'  Prep {i+1}/{n_total} ({len(samples)} quads)', flush=True)
    print(f'  Quads found: {len(samples)}/{n_total}')
    return samples


def _process_quad(args):
    idx, sample, cfg = args
    theta, theta_E, bx, by, macro_result = sample
    rng = np.random.default_rng(cfg['seed_base'] + idx)
    macro = MacroLens()
    _ = macro.build(theta['z_l'], theta['z_s'],
        theta['sigma_v'], theta['q'], theta['phi_lens'],
        theta['gamma_ext'], theta['theta_gamma'])
    if cfg.get('los'):
        los_pop = LOSPopulation()
        los_kwargs = los_pop.realise(theta['z_l'], theta['z_s'], rng)
        macro.add_los(los_kwargs)
    if cfg.get('substructure'):
        sub_pop = SubhaloPopulation(f_sub=cfg.get('f_sub'))
        subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
        macro.add_substructure(subhalos)
        full_result = macro.solve(bx, by, num_random=2, search_window=2)
        if full_result is None or full_result['n_images'] < 4:
            return None
        theta_x, theta_y, mu = full_result['theta_x'], full_result['theta_y'], full_result['mu']
    else:
        theta_x, theta_y, mu = macro_result['theta_x'], macro_result['theta_y'], macro_result['mu']
    if cfg.get('selection'):
        if not passes_selection(theta_x, theta_y, mu, theta['source_flux']):
            return None
    F_true = mu * theta['source_flux']
    if cfg.get('noise'):
        F_obs = RadioNoise().apply(np.abs(F_true), rng)
    else:
        F_obs = np.abs(F_true)
    r = compute_rfold(theta_x, theta_y, mu, F_obs)
    return r if r >= 0 else None


def run_config(samples, config, n_max=None, n_workers=1):
    if n_max is not None:
        samples = samples[:n_max]
    n_total = len(samples)
    args_list = [(i, s, config) for i, s in enumerate(samples)]
    if n_workers <= 1:
        R_fold = []
        for idx, sample in enumerate(samples):
            r = _process_quad((idx, sample, config))
            if r is not None:
                R_fold.append(r)
            if (idx + 1) % 100 == 0:
                print(f'    quad {idx+1}/{n_total} ({len(R_fold)} valid)', flush=True)
        return np.array(R_fold)
    batches = [args_list[i:i + n_workers * 2] for i in range(0, len(args_list), n_workers * 2)]
    R_fold = []
    done = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for batch in batches:
            for r in pool.map(_process_quad, batch):
                if r is not None:
                    R_fold.append(r)
            done += len(batch)
            print(f'    {done}/{n_total} ({len(R_fold)} valid)', flush=True)
    return np.array(R_fold)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-systems', type=int, default=500)
    parser.add_argument('--parallel', action='store_true')
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--fast', action='store_true', help='Alias for --n-systems 200 --parallel')
    args = parser.parse_args()

    if args.fast:
        args.n_systems = 200
        args.parallel = True
    n_workers = args.workers or os.cpu_count() if args.parallel else 1

    outputs_dir = Path(__file__).parent / 'outputs'
    outputs_dir.mkdir(exist_ok=True)

    from data.radio_quads import build_observed_catalog as get_obs
    R_obs = []
    for s in get_obs():
        r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
            np.array(s['parity']), np.array(s['fluxes_mjy']))
        if r >= 0: R_obs.append(r)
    R_obs = np.array(R_obs)
    print(f'Observed: {len(R_obs)} systems')

    print('Preparing lenses...')
    samples = prepare_lenses(args.n_systems, 42)

    f_sub_values = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02]
    base_cfg = {'substructure': True, 'los': True, 'selection': True, 'noise': True}
    results = {}

    print('=== f_sub grid ===')
    for fs in f_sub_values:
        ckpt = outputs_dir / f'fsub_{fs:.3f}.npy'
        if ckpt.exists():
            R_sim = np.load(ckpt)
            print(f'  f_sub={fs:.3f} (cached, N={len(R_sim)})')
        else:
            cfg_run = {**base_cfg, 'f_sub': fs, 'seed_base': 42 + int(fs * 1000)}
            t0 = time.perf_counter()
            R_sim = run_config(samples, cfg_run, n_workers=n_workers)
            np.save(ckpt, R_sim)
            print(f'  f_sub={fs:.3f}: N={len(R_sim)}, {time.perf_counter()-t0:.0f}s')
        results[fs] = R_sim
        if len(R_sim) > 0:
            res = compare(R_obs, R_sim, n_bootstrap=1000)
            print(f'    KS_p={res["ks_p_value"]:.4f}  Wass={res["wasserstein_distance"]:.4f}  Mean={np.mean(R_sim):.4f}')

    print('\n=== Posterior ===')
    from src.inference import run_inference
    f_grid = [fs for fs in f_sub_values if len(results[fs]) > 0]
    R_grid = [results[fs] for fs in f_grid]
    post = run_inference(f_grid, R_grid, R_obs, prior='uniform')
    print(f'  Best-fit f_sub = {post["best_fit"]:.4f}')
    print(f'  Median   f_sub = {post["median"]:.4f}')
    print(f'  68% CI        = [{post["ci_68_low"]:.4f}, {post["ci_68_high"]:.4f}]')
    print(f'  Mean +- std   = {post["mean"]:.4f} +- {post["std"]:.4f}')
    np.savez(outputs_dir / 'fsub_posterior.npz',
        f_sub_vals=post['f_sub_vals'], posterior=post['posterior'],
        best_fit=np.array([post['best_fit']]), median=np.array([post['median']]),
        ci_68_low=np.array([post['ci_68_low']]), ci_68_high=np.array([post['ci_68_high']]))


if __name__ == '__main__':
    main()
