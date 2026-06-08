import sys
import numpy as np
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation, LOSPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from src.statistic import compute_rfold, compute_cusp_statistic
from src.comparison import compare
from data.radio_quads import build_observed_catalog

RNG_SEED = 42
CONFIGS = [
    {"name": "macro_only",         "substructure": False, "los": False, "selection": False, "noise": False},
    {"name": "macro_los",          "substructure": False, "los": True,  "selection": False, "noise": False},
    {"name": "macro_subhalos",     "substructure": True,  "los": False, "selection": False, "noise": False},
    {"name": "macro_sub_los",      "substructure": True,  "los": True,  "selection": False, "noise": False},
    {"name": "macro_sub_sel",      "substructure": True,  "los": False, "selection": True,  "noise": False},
    {"name": "full",               "substructure": True,  "los": False, "selection": True,  "noise": True},
]


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


def prepare_lenses(n_total, seed):
    rng = np.random.default_rng(seed)
    pop = LensPopulation(rng)
    samples = []
    for i in range(n_total):
        theta = pop.sample(rng)
        macro = MacroLens()
        theta_E = macro.build(
            theta['z_l'], theta['z_s'],
            theta['sigma_v'], theta['q'], theta['phi_lens'],
            theta['gamma_ext'], theta['theta_gamma'],
        )
        bx, by, macro_result = sample_source_in_caustic(
            rng, theta_E, theta['q'], macro)
        if macro_result is not None:
            samples.append((theta, theta_E, bx, by, macro_result))
        if (i + 1) % 500 == 0:
            print(f'  Prep {i+1}/{n_total} ({len(samples)} quads)', flush=True)
    print(f'  Quads found: {len(samples)}/{n_total}')
    return samples


def run_config(samples, config, rng, name=''):
    f_sub = config.get('f_sub', None)
    sub_pop = SubhaloPopulation(f_sub=f_sub)
    los_pop = LOSPopulation()
    radio_noise = RadioNoise()

    R_fold, R_cusp = [], []
    n_total = len(samples)

    for idx, (theta, theta_E, bx, by, macro_result) in enumerate(samples):
        if (idx + 1) % 200 == 0:
            print(f'    {name}: {idx+1}/{n_total} ({len(R_fold)} valid)', flush=True)
        macro = MacroLens()
        _ = macro.build(
            theta['z_l'], theta['z_s'],
            theta['sigma_v'], theta['q'], theta['phi_lens'],
            theta['gamma_ext'], theta['theta_gamma'],
        )

        if config['los']:
            _ = los_pop.realise(theta['z_l'], theta['z_s'], rng)

        if config['substructure']:
            subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
            macro.add_substructure(subhalos)
            full_result = macro.solve(bx, by, num_random=2, search_window=2)
            if full_result is None or full_result['n_images'] < 4:
                continue
            theta_x = full_result['theta_x']
            theta_y = full_result['theta_y']
            mu = full_result['mu']
        else:
            theta_x = macro_result['theta_x']
            theta_y = macro_result['theta_y']
            mu = macro_result['mu']

        if config['selection']:
            if not passes_selection(theta_x, theta_y, mu,
                                    theta['source_flux']):
                continue

        F_true = mu * theta['source_flux']
        if config['noise']:
            F_obs = radio_noise.apply(np.abs(F_true), rng)
        else:
            F_obs = np.abs(F_true)

        r = compute_rfold(theta_x, theta_y, mu, F_obs)
        if r >= 0:
            R_fold.append(r)
        # cusp: use signed mu directly (no F_obs) for simulated data
        c = compute_cusp_statistic(theta_x, theta_y, mu)
        if c >= 0:
            R_cusp.append(c)

    return np.array(R_fold), np.array(R_cusp)


def load_observed():
    systems = build_observed_catalog()
    R_fold, R_cusp = [], []
    for s in systems:
        t_x = np.array(s['x_arcsec'])
        t_y = np.array(s['y_arcsec'])
        p = np.array(s['parity'])
        f = np.array(s['fluxes_mjy'])
        r = compute_rfold(t_x, t_y, p, f)
        if r >= 0:
            R_fold.append(r)
        c = compute_cusp_statistic(t_x, t_y, p, f)
        if c >= 0:
            R_cusp.append(c)
    return np.array(R_fold), np.array(R_cusp), systems


def run_all_ablations(n_total, outputs_dir):
    print('=== Preparing shared lens sample ===')
    samples = prepare_lenses(n_total, RNG_SEED)
    R_fold_obs, R_cusp_obs, systems = load_observed()

    print(f'\nObserved: N_Rfold={len(R_fold_obs)}, N_Rcusp={len(R_cusp_obs)}')
    print(f'R_fold_obs = {np.round(R_fold_obs, 3)}')
    print(f'R_cusp_obs = {np.round(R_cusp_obs, 3)}')
    print()

    print('=== Ablation matrix ===')
    rows = []
    for cfg in CONFIGS:
        name = cfg['name']
        print(f'  Config: {name}', end='', flush=True)
        ckpt_path = outputs_dir / f'ablation_{name}.npz'
        if ckpt_path.exists():
            data = np.load(ckpt_path)
            R_sim = data['R_fold']
            R_cusp_sim = data['R_cusp']
            print(f' (cached, N={len(R_sim)})')
        else:
            rng = np.random.default_rng(RNG_SEED + hash(name) % 2**31)
            R_sim, R_cusp_sim = run_config(samples, cfg, rng, name=name)
            np.savez(ckpt_path, R_fold=R_sim, R_cusp=R_cusp_sim)
            print(f' (N={len(R_sim)})')

        if len(R_sim) > 0:
            result = compare(R_fold_obs, R_sim, n_bootstrap=1000)
            rows.append({
                'name': name,
                'N_sim': len(R_sim),
                'ks_stat': result['ks_statistic'],
                'ks_p': result['ks_p_value'],
                'wass': result['wasserstein_distance'],
                'mean_r': float(np.mean(R_sim)),
                'var_r': float(np.var(R_sim)),
            })
        else:
            rows.append({
                'name': name, 'N_sim': 0,
                'ks_stat': -1, 'ks_p': -1, 'wass': -1,
                'mean_r': -1, 'var_r': -1,
            })

    print()
    print('+' + '-' * 97 + '+')
    header = '| {:<17} | {:>6} | {:>9} | {:>9} | {:>9} | {:>7} | {:>9} |'.format(
        'Model', 'N_sim', 'KS_stat', 'KS_p', 'Wass', 'Mean_R', 'Var_R')
    print(header)
    print('|' + '-' * 97 + '|')
    for r in rows:
        if r['N_sim'] > 0:
            line = '| {:<17} | {:>6d} | {:>9.4f} | {:>9.4f} | {:>9.4f} | {:>7.4f} | {:>9.6f} |'.format(
                r['name'], r['N_sim'], r['ks_stat'], r['ks_p'],
                r['wass'], r['mean_r'], r['var_r'])
        else:
            line = '| {:<17} | {:>6d} | {:>9} | {:>9} | {:>9} | {:>7} | {:>9} |'.format(
                r['name'], 0, '-', '-', '-', '-', '-')
        print(line)
    print('+' + '-' * 97 + '+')
    print()

    print('=== Cusp statistic ===')
    cusp_rows = []
    for cfg in CONFIGS:
        name = cfg['name']
        ckpt_path = outputs_dir / f'ablation_{name}.npz'
        if ckpt_path.exists():
            data = np.load(ckpt_path)
            R_cusp_sim = data['R_cusp']
        else:
            continue
        if len(R_cusp_sim) > 0 and len(R_cusp_obs) > 0:
            result = compare(R_cusp_obs, R_cusp_sim, n_bootstrap=1000)
            cusp_rows.append({
                'name': name,
                'N_sim': len(R_cusp_sim),
                'ks_stat': result['ks_statistic'],
                'ks_p': result['ks_p_value'],
                'wass': result['wasserstein_distance'],
                'mean_c': float(np.mean(R_cusp_sim)),
            })

    if cusp_rows:
        print('+' + '-' * 89 + '+')
        header = '| {:<17} | {:>6} | {:>9} | {:>9} | {:>9} | {:>9} |'.format(
            'Model', 'N_sim', 'KS_stat', 'KS_p', 'Wass', 'Mean_Rc')
        print(header)
        print('|' + '-' * 89 + '|')
        for r in cusp_rows:
            line = '| {:<17} | {:>6d} | {:>9.4f} | {:>9.4f} | {:>9.4f} | {:>9.6f} |'.format(
                r['name'], r['N_sim'], r['ks_stat'], r['ks_p'],
                r['wass'], r['mean_c'])
            print(line)
        print('+' + '-' * 89 + '+')
    print()

    return rows, samples


def leave_one_out(R_sim, systems, R_cusp_sim=None):
    print('=== Leave-one-out sensitivity ===')
    print(f'{"Left out":<14} {"R_fold":>8} {"KS_p":>9} {"dKS_p":>9}')
    print('-' * 44)
    full_ks_p = compare(
        np.array([compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                                np.array(s['parity']), np.array(s['fluxes_mjy']))
                  for s in systems if compute_rfold(np.array(s['x_arcsec']),
                                                    np.array(s['y_arcsec']),
                                                    np.array(s['parity']),
                                                    np.array(s['fluxes_mjy'])) >= 0]),
        R_sim)['ks_p_value']

    for exclude in systems:
        n_excl = exclude['name']
        included = [s for s in systems if s['name'] != n_excl]
        R_loo = []
        for s in included:
            r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                              np.array(s['parity']), np.array(s['fluxes_mjy']))
            if r >= 0:
                R_loo.append(r)
        if len(R_loo) > 0:
            result = compare(np.array(R_loo), R_sim, n_bootstrap=1000)
            delta_p = result['ks_p_value'] - full_ks_p
            r_val = compute_rfold(np.array(exclude['x_arcsec']),
                                  np.array(exclude['y_arcsec']),
                                  np.array(exclude['parity']),
                                  np.array(exclude['fluxes_mjy']))
            print(f'{n_excl:<14} {r_val:>8.3f} {result["ks_p_value"]:>9.4f} {delta_p:>+9.4f}')
    print()


def parity_perturbation(R_sim, systems):
    print('=== Parity perturbation test ===')
    print('Flipping one parity element at a time:')
    print(f'{"System":<14} {"Flip":>6} {"R_fold":>8} {"KS_p":>9} {"dKS_p":>9}')
    print('-' * 50)

    R_obs_all = np.array([
        compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                      np.array(s['parity']), np.array(s['fluxes_mjy']))
        for s in systems
    ])
    # Filter valid
    valid_mask = R_obs_all >= 0
    full_ks_p = compare(R_obs_all[valid_mask], R_sim)['ks_p_value']

    for s in systems:
        base_parity = np.array(s['parity'])
        for flip_idx in range(4):
            perturbed = base_parity.copy()
            perturbed[flip_idx] *= -1
            r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                              perturbed, np.array(s['fluxes_mjy']))
            if r < 0:
                continue
            # Recompute with perturbed parity for this system
            R_pert = []
            for s2 in systems:
                p2 = np.array(s2['parity'])
                if s2['name'] == s['name']:
                    p2 = perturbed
                r2 = compute_rfold(np.array(s2['x_arcsec']),
                                   np.array(s2['y_arcsec']),
                                   p2, np.array(s2['fluxes_mjy']))
                if r2 >= 0:
                    R_pert.append(r2)
            if len(R_pert) > 0:
                result = compare(np.array(R_pert), R_sim, n_bootstrap=1000)
                delta_p = result['ks_p_value'] - full_ks_p
                print(f'{s["name"]:<14} {flip_idx:>6d} {r:>8.3f} '
                      f'{result["ks_p_value"]:>9.4f} {delta_p:>+9.4f}')
    print()

    print('Monte Carlo: random single-element flip per system, 1000 iterations')
    rng = np.random.default_rng(999)
    p_dist = []
    for _ in range(1000):
        R_pert = []
        for s in systems:
            p = np.array(s['parity'])
            flip_i = rng.integers(0, 4)
            p[flip_i] *= -1
            r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                              p, np.array(s['fluxes_mjy']))
            if r >= 0:
                R_pert.append(r)
        if len(R_pert) > 0:
            res = compare(np.array(R_pert), R_sim, n_bootstrap=1000)
            p_dist.append(res['ks_p_value'])
    p_dist = np.array(p_dist)
    print(f'  KS_p under perturbation: median={np.median(p_dist):.4f}, '
          f'p<0.05={np.mean(p_dist < 0.05)*100:.0f}%, '
          f'p<0.01={np.mean(p_dist < 0.01)*100:.0f}%')
    print()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-systems', type=int, default=None,
                        help='Override n_systems from config (for faster ablation)')
    args = parser.parse_args()

    config_path = Path(__file__).parent.parent / 'config.yaml'
    with open(config_path, encoding='utf-8') as f:
        cfg = yaml.safe_load(f)

    n_total = args.n_systems or cfg['simulation']['n_systems']
    outputs_dir = Path(__file__).parent.parent / 'outputs'
    outputs_dir.mkdir(exist_ok=True)

    # Step 1: run all ablation configurations
    rows, samples = run_all_ablations(n_total, outputs_dir)

    # Step 2: leave-one-out against full model R_sim
    ckpt_full = outputs_dir / 'ablation_full.npz'
    if ckpt_full.exists():
        R_sim_full = np.load(ckpt_full)['R_fold']
    else:
        print('Full model results not found — re-run main pipeline first')
        return

    _, _, systems = load_observed()
    leave_one_out(R_sim_full, systems)
    parity_perturbation(R_sim_full, systems)

    # Step 3: f_sub response surface
    print('=== f_sub response surface ===')
    f_sub_values = [0.0, 0.002, 0.005, 0.01, 0.02, 0.05]
    print(f'{"f_sub":>8} {"KS_p":>9} {"Wass":>9} {"Mean_R":>8}')
    print('-' * 38)

    f_sub_rng = np.random.default_rng(RNG_SEED + 9999)
    for fs in f_sub_values:
        ckpt_path = outputs_dir / f'ablation_fsub_{fs:.3f}.npz'
        if ckpt_path.exists():
            data = np.load(ckpt_path)
            R_sim_fs = data['R_fold']
        else:
            R_sim_fs, _ = run_config(samples, {
                "substructure": True, "los": False,
                "selection": True, "noise": True,
                "f_sub": fs,
            }, np.random.default_rng(f_sub_rng.integers(0, 2**31)),
                name=f'f_sub={fs:.3f}')

            np.savez(ckpt_path, R_fold=R_sim_fs)

        if len(R_sim_fs) > 0 and len(R_fold_obs := load_observed()[0]) > 0:
            res = compare(R_fold_obs, R_sim_fs, n_bootstrap=1000)
            print(f'{fs:>8.3f} {res["ks_p_value"]:>9.4f} '
                  f'{res["wasserstein_distance"]:>9.4f} '
                  f'{np.mean(R_sim_fs):>8.4f}')
    print()


if __name__ == '__main__':
    main()
