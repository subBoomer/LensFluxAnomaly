import sys
import numpy as np
import yaml
from pathlib import Path

from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from src.statistic import compute_rfold
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


def load_checkpoint(outputs_dir):
    ckpt_path = outputs_dir / 'checkpoint.npz'
    if ckpt_path.exists():
        data = np.load(ckpt_path)
        return int(data['completed'].item()), data['R_sim'].tolist()
    return 0, []


def save_checkpoint(outputs_dir, completed, R_sim):
    outputs_dir.mkdir(exist_ok=True)
    np.savez(outputs_dir / 'checkpoint.npz',
             completed=np.array([completed]),
             R_sim=np.array(R_sim))


def main():
    config_path = Path(__file__).parent / 'config.yaml'
    with open(config_path, encoding='utf-8') as f:
        config = yaml.safe_load(f)

    observed_systems = get_observed_catalog()
    R_obs = []
    for sys in observed_systems:
        parity = np.array(sys.get('parity', [1, 1, 1, 1]))
        r = compute_rfold(
            np.array(sys['x_arcsec']),
            np.array(sys['y_arcsec']),
            parity,
            np.array(sys['fluxes_mjy']),
        )
        if r >= 0:
            R_obs.append(r)
    R_obs = np.array(R_obs)

    cfg = config['simulation']
    n_total = cfg['n_systems']
    ckpt_interval = cfg.get('checkpoint_interval', 100)

    outputs_dir = Path(__file__).parent / 'outputs'
    completed, R_sim = load_checkpoint(outputs_dir)
    if completed >= n_total:
        print('Inference already completed. Use --force to re-run.')
        return

    rng = np.random.default_rng(cfg['seed'])
    pop = LensPopulation(rng)
    radio_noise = RadioNoise()
    sub_pop = SubhaloPopulation()

    for _ in range(completed):
        theta = pop.sample(rng)
        _ = theta

    n_quad = 0
    n_selected = 0
    start = completed

    for i in range(start, n_total):
        theta = pop.sample(rng)
        macro = MacroLens()
        theta_E = macro.build(
            theta['z_l'], theta['z_s'],
            theta['sigma_v'], theta['q'], theta['phi_lens'],
            theta['gamma_ext'], theta['theta_gamma'],
        )
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

        if (i - start + 1) % ckpt_interval == 0:
            save_checkpoint(outputs_dir, i + 1, R_sim)
            elapsed = (i - start + 1) / n_total * 100
            print(f'  {i+1}/{n_total} ({elapsed:.0f}%): quad={n_quad}, selected={n_selected}, R_fold_values={len(R_sim)}', flush=True)

    R_sim = np.array(R_sim)
    result = compare(R_obs, R_sim, config['comparison']['n_bootstrap'])

    save_dict = {
        'R_obs': R_obs,
        'R_sim': R_sim,
        'n_observed': np.array([result['n_observed']]),
        'n_simulated': np.array([result['n_simulated']]),
        'ks_statistic': np.array([result['ks_statistic']]),
        'ks_p_value': np.array([result['ks_p_value']]),
        'wasserstein_distance': np.array([result['wasserstein_distance']]),
    }
    np.savez(outputs_dir / 'inference_result.npz', **save_dict)
    Path(outputs_dir / 'checkpoint.npz').unlink(missing_ok=True)

    print(f'Observed systems: {len(R_obs)}')
    print(f'Simulated systems: {len(R_sim)}')
    print(f'KS statistic: {result["ks_statistic"]:.4f} (p={result["ks_p_value"]:.4f})')
    print(f'Wasserstein distance: {result["wasserstein_distance"]:.4f}')


if __name__ == '__main__':
    main()
