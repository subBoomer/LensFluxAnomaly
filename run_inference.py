import argparse
import os
import sys
import time
import numpy as np
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from src.config import validate as validate_config
from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation, LOSPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from src.statistic import compute_rfold, compute_cusp_statistic
from src.comparison import compare
from data.radio_quads import build_observed_catalog as get_observed_catalog


def sample_source_in_caustic(rng, theta_E, q, macro, max_attempts=30):
    caustic_size = theta_E * (1.0 - q) / (1.0 + q)
    scale = max(caustic_size * 0.8, 0.01)
    for _ in range(max_attempts):
        bx = rng.uniform(-scale, scale)
        by = rng.uniform(-scale, scale)
        result = macro.solve_macro_only(bx, by)
        if result is not None and result['n_images'] == 4:
            return bx, by, result
    return None, None, None


def simulate_one(seed, params):
    rng = np.random.default_rng(seed)
    pop = LensPopulation(rng)
    sub_pop = SubhaloPopulation(f_sub=params.get('f_sub'))
    los_pop = LOSPopulation()
    radio_noise = RadioNoise()
    radio_noise.snr = params['radio_snr']

    theta = pop.sample(rng)
    macro = MacroLens()
    theta_E = macro.build(
        theta['z_l'], theta['z_s'],
        theta['sigma_v'], theta['q'], theta['phi_lens'],
        theta['gamma_ext'], theta['theta_gamma'],
    )
    los_kwargs = los_pop.realise(theta['z_l'], theta['z_s'], rng)
    macro.add_los(los_kwargs)
    bx, by, macro_result = sample_source_in_caustic(rng, theta_E, theta['q'], macro)
    if macro_result is None:
        return None, None, 0, 0
    if not passes_selection(
        macro_result['theta_x'], macro_result['theta_y'],
        macro_result['mu'], theta['source_flux'],
    ):
        return None, None, 1, 0
    subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
    macro.add_substructure(subhalos)
    full_result = macro.solve(bx, by, num_random=2, search_window=2)
    if full_result is None or full_result['n_images'] < 4:
        return None, None, 1, 1
    F_true = full_result['mu'] * theta['source_flux']
    F_obs = radio_noise.apply(F_true, rng)
    R_fold = compute_rfold(
        full_result['theta_x'], full_result['theta_y'],
        full_result['mu'], F_obs,
    )
    if R_fold < 0:
        return None, None, 1, 1
    R_cusp = compute_cusp_statistic(
        full_result['theta_x'], full_result['theta_y'],
        full_result['mu'], None,
    )
    return R_fold, R_cusp, 1, 1


def simulate_batch(seeds, params):
    fold_vals = []
    cusp_vals = []
    n_quad = 0
    n_sel = 0
    for seed in seeds:
        r, c, q, s = simulate_one(seed, params)
        n_quad += q
        n_sel += s
        if r is not None:
            fold_vals.append(r)
            cusp_vals.append(c)
    return fold_vals, cusp_vals, n_quad, n_sel


def load_checkpoint(outputs_dir):
    ckpt_path = outputs_dir / 'checkpoint.npz'
    if ckpt_path.exists():
        data = np.load(ckpt_path)
        R_sim = data['R_sim'].tolist() if 'R_sim' in data else data['R_fold'].tolist()
        C_sim = data['C_sim'].tolist() if 'C_sim' in data else []
        return int(data['completed'].item()), R_sim, C_sim
    return 0, [], []


def save_checkpoint(outputs_dir, completed, R_sim, C_sim):
    outputs_dir.mkdir(exist_ok=True)
    np.savez(outputs_dir / 'checkpoint.npz',
             completed=np.array([completed]),
             R_sim=np.array(R_sim),
             C_sim=np.array(C_sim))


def main():
    parser = argparse.ArgumentParser(description='LensFluxAnomaly inference pipeline')
    parser.add_argument('--force', action='store_true', help='Re-run even if results exist')
    parser.add_argument('--n-systems', type=int, default=None, help='Override n_systems from config')
    parser.add_argument('--parallel', action='store_true', help='Use multiprocessing')
    parser.add_argument('--workers', type=int, default=None, help='Number of worker processes (default: CPU count)')
    parser.add_argument('--quick', action='store_true', help='Fast dev mode: 50 draws, force re-run')
    parser.add_argument('--snr', type=float, default=None, help='Override radio_snr from config')
    args = parser.parse_args()

    config_path = Path(__file__).parent / 'config.yaml'
    with open(config_path, encoding='utf-8') as f:
        config = yaml.safe_load(f)
    validate_config(config)

    if args.n_systems is not None:
        config['simulation']['n_systems'] = args.n_systems
    if args.quick:
        config['simulation']['n_systems'] = 50
        args.force = True
        args.parallel = True
    if args.snr is not None:
        config['noise']['radio_snr'] = args.snr

    observed_systems = get_observed_catalog()
    R_obs, C_obs = [], []
    for sys in observed_systems:
        parity = np.array(sys.get('parity', [1, 1, 1, 1]))
        fluxes = np.array(sys['fluxes_mjy'])
        x = np.array(sys['x_arcsec'])
        y = np.array(sys['y_arcsec'])
        r = compute_rfold(x, y, parity, fluxes)
        if r >= 0:
            R_obs.append(r)
            c = compute_cusp_statistic(x, y, parity, fluxes)
            if c >= 0:
                C_obs.append(c)
    R_obs = np.array(R_obs)
    C_obs = np.array(C_obs)

    cfg = config['simulation']
    n_total = cfg['n_systems']
    ckpt_interval = cfg.get('checkpoint_interval', 100)

    outputs_dir = Path(__file__).parent / 'outputs'
    completed, R_sim, C_sim = load_checkpoint(outputs_dir)
    if args.force:
        if completed > 0 or (outputs_dir / 'inference_result.npz').exists():
            print('Forcing re-run (clearing previous results)')
        completed = 0
        R_sim, C_sim = [], []
    elif completed >= n_total:
        print('Inference already completed. Use --force to re-run.')
        return

    params = {
        'radio_snr': config['noise']['radio_snr'],
        'f_sub': config['substructure']['f_sub_dmo'],
    }

    n_quad = 0
    n_selected = 0
    start = completed

    if args.parallel:
        n_workers = args.workers or os.cpu_count()
        batch_size = max(1, ckpt_interval)
        seeds = list(range(start, n_total))
        batches = [seeds[i:i + batch_size] for i in range(0, len(seeds), batch_size)]

        print(f'Running {n_total - start} draws on {n_workers} workers ({len(batches)} batches of ~{batch_size})')
        t_start = time.perf_counter()

        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            fut_to_batch = {pool.submit(simulate_batch, b, params): b for b in batches}
            done = 0
            for fut in as_completed(fut_to_batch):
                batch_fold, batch_cusp, b_quad, b_sel = fut.result()
                R_sim.extend(batch_fold)
                C_sim.extend(batch_cusp)
                n_quad += b_quad
                n_selected += b_sel
                done += len(fut_to_batch[fut])
                completed = start + done
                if completed % ckpt_interval == 0 or done == n_total - start:
                    save_checkpoint(outputs_dir, completed, R_sim, C_sim)
                    elapsed = time.perf_counter() - t_start
                    rate = done / elapsed if elapsed > 0 else 0
                    print(f'  {completed}/{n_total} ({done/elapsed:.0f} draws/s): quad={n_quad}, selected={n_selected}, R_fold={len(R_sim)}, R_cusp={len(C_sim)}', flush=True)

    else:
        rng = np.random.default_rng(cfg['seed'])
        pop = LensPopulation(rng)
        radio_noise = RadioNoise()
        sub_pop = SubhaloPopulation()

        for _ in range(completed):
            theta = pop.sample(rng)

        los_pop = LOSPopulation()
        for i in range(start, n_total):
            theta = pop.sample(rng)
            macro = MacroLens()
            theta_E = macro.build(
                theta['z_l'], theta['z_s'],
                theta['sigma_v'], theta['q'], theta['phi_lens'],
                theta['gamma_ext'], theta['theta_gamma'],
            )
            los_kwargs = los_pop.realise(theta['z_l'], theta['z_s'], rng)
            macro.add_los(los_kwargs)
            bx, by, macro_result = sample_source_in_caustic(rng, theta_E, theta['q'], macro)
            if macro_result is None:
                continue
            n_quad += 1
            if not passes_selection(
                macro_result['theta_x'], macro_result['theta_y'],
                macro_result['mu'], theta['source_flux'],
            ):
                continue
            n_selected += 1
            subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
            macro.add_substructure(subhalos)
            full_result = macro.solve(bx, by, num_random=2, search_window=2)
            if full_result is None or full_result['n_images'] < 4:
                continue
            F_true = full_result['mu'] * theta['source_flux']
            F_obs = radio_noise.apply(F_true, rng)
            R_fold = compute_rfold(
                full_result['theta_x'], full_result['theta_y'],
                full_result['mu'], F_obs,
            )
            if R_fold >= 0:
                R_sim.append(R_fold)
                R_cusp = compute_cusp_statistic(
                    full_result['theta_x'], full_result['theta_y'],
                    full_result['mu'], None,
                )
                if R_cusp >= 0:
                    C_sim.append(R_cusp)

            if (i - start + 1) % ckpt_interval == 0:
                save_checkpoint(outputs_dir, i + 1, R_sim, C_sim)
                elapsed = (i - start + 1) / n_total * 100
                print(f'  {i+1}/{n_total} ({elapsed:.0f}%): quad={n_quad}, selected={n_selected}, R_fold={len(R_sim)}, R_cusp={len(C_sim)}', flush=True)

    R_sim = np.array(R_sim)
    C_sim = np.array(C_sim)
    result = compare(R_obs, R_sim, config['comparison']['n_bootstrap'])
    cusp_result = compare(C_obs, C_sim, config['comparison']['n_bootstrap']) if len(C_sim) >= 5 and len(C_obs) >= 5 else None

    save_dict = {
        'R_obs': R_obs,
        'R_sim': R_sim,
        'n_observed': np.array([result['n_observed']]),
        'n_simulated': np.array([result['n_simulated']]),
        'ks_statistic': np.array([result['ks_statistic']]),
        'ks_p_value': np.array([result['ks_p_value']]),
        'wasserstein_distance': np.array([result['wasserstein_distance']]),
        'cusp_obs': C_obs,
        'cusp_sim': C_sim,
    }
    if cusp_result is not None:
        save_dict['cusp_ks_statistic'] = np.array([cusp_result['ks_statistic']])
        save_dict['cusp_ks_p_value'] = np.array([cusp_result['ks_p_value']])
        save_dict['cusp_wasserstein_distance'] = np.array([cusp_result['wasserstein_distance']])
    np.savez(outputs_dir / 'inference_result.npz', **save_dict)
    Path(outputs_dir / 'checkpoint.npz').unlink(missing_ok=True)

    print(f'Observed systems: {len(R_obs)}')
    print(f'Simulated systems: {len(R_sim)}')
    print(f'R_fold — KS statistic: {result["ks_statistic"]:.4f} (p={result["ks_p_value"]:.4f})')
    print(f'R_fold — Wasserstein distance: {result["wasserstein_distance"]:.4f}')
    if cusp_result is not None:
        print(f'R_cusp — KS statistic: {cusp_result["ks_statistic"]:.4f} (p={cusp_result["ks_p_value"]:.4f})')
        print(f'R_cusp — Wasserstein distance: {cusp_result["wasserstein_distance"]:.4f}')


if __name__ == '__main__':
    main()
