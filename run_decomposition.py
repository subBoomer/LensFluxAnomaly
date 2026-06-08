import sys, numpy as np, time, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation, LOSPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from src.statistic import compute_rfold, compute_cusp_statistic


def prepare_macro_quads(n_total, seed):
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
        cs = theta_E * (1.0 - theta['q']) / (1.0 + theta['q'])
        scale = max(cs * 0.8, 0.01)
        for _ in range(30):
            bx = rng.uniform(-scale, scale)
            by = rng.uniform(-scale, scale)
            result = macro.solve_macro_only(bx, by)
            if result is not None and result['n_images'] == 4:
                samples.append((theta, theta_E, bx, by, result))
                break
        if (i + 1) % 200 == 0:
            print(f'  Prep {i+1}/{n_total} ({len(samples)} quads)', flush=True)
    print(f'  Quads found: {len(samples)}/{n_total}')
    return samples


def _solve_and_compute(macro, bx, by, source_flux, snr, rng):
    full_result = macro.solve(bx, by, num_random=2, search_window=2)
    if full_result is None or full_result['n_images'] < 4:
        return None, None
    if not passes_selection(
        full_result['theta_x'], full_result['theta_y'],
        full_result['mu'], source_flux, rng,
    ):
        return None, None
    F_true = full_result['mu'] * source_flux
    F_obs = RadioNoise().apply(np.abs(F_true), rng)
    R_fold = compute_rfold(
        full_result['theta_x'], full_result['theta_y'],
        full_result['mu'], F_obs,
    )
    if R_fold < 0:
        return None, None
    R_cusp = compute_cusp_statistic(
        full_result['theta_x'], full_result['theta_y'],
        full_result['mu'], None,
    )
    return R_fold, R_cusp


def _process_cell(args):
    cell_idx, samples, config = args
    rng = np.random.default_rng(config['seed'] + cell_idx)
    los_pop = LOSPopulation()
    sub_pop = SubhaloPopulation(f_sub=config['f_sub'])
    snr = config.get('snr', 20)

    folds, cusps = [], []
    for idx, (theta, theta_E, bx, by, macro_result) in enumerate(samples):
        macro = MacroLens()
        macro.build(
            theta['z_l'], theta['z_s'],
            theta['sigma_v'], theta['q'], theta['phi_lens'],
            theta['gamma_ext'], theta['theta_gamma'],
        )
        if config['los']:
            los_kwargs = los_pop.realise(theta['z_l'], theta['z_s'], rng)
            macro.add_los(los_kwargs)
        if config['substructure']:
            subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
            macro.add_substructure(subhalos)
            r_fold, r_cusp = _solve_and_compute(macro, bx, by, theta['source_flux'], snr, rng)
        else:
            theta_x = macro_result['theta_x']
            theta_y = macro_result['theta_y']
            mu = macro_result['mu']
            if not passes_selection(theta_x, theta_y, mu, theta['source_flux'], rng):
                continue
            F_true = mu * theta['source_flux']
            F_obs = RadioNoise().apply(np.abs(F_true), rng)
            r_fold = compute_rfold(theta_x, theta_y, mu, F_obs)
            if r_fold < 0:
                continue
            r_cusp = compute_cusp_statistic(theta_x, theta_y, mu, None)
        if r_fold is not None:
            folds.append(r_fold)
            cusps.append(r_cusp)
    return folds, cusps


def run_decomposition_table(quads, f_sub_list, configs, n_workers):
    rows = []
    for cfg_name, cfg in configs.items():
        for fs in f_sub_list:
            cell_config = {**cfg, 'f_sub': fs, 'seed': 42}
            t0 = time.perf_counter()
            if cfg['substructure'] and n_workers > 1:
                chunk_size = max(1, len(quads) // n_workers)
                chunks = [quads[i:i + chunk_size] for i in range(0, len(quads), chunk_size)]
                args_list = [(i, ch, cell_config) for i, ch in enumerate(chunks)]
                all_folds, all_cusps = [], []
                with ProcessPoolExecutor(max_workers=n_workers) as pool:
                    for folds, cusps in pool.map(_process_cell, args_list):
                        all_folds.extend(folds)
                        all_cusps.extend(cusps)
                fold_arr = np.array(all_folds)
                cusp_arr = np.array(all_cusps)
            else:
                fold_arr, cusp_arr = _process_cell((0, quads, cell_config))
                fold_arr = np.array(fold_arr)
                cusp_arr = np.array(cusp_arr)
            dt = time.perf_counter() - t0
            label = f'{cfg_name}  f_sub={fs}'
            if len(fold_arr) < 5:
                rows.append((label, 0, 0, 0, 0, 0, 0, dt))
                print(f'  {label}: N={len(fold_arr)} <- too few, skipped', flush=True)
                continue
            n = len(fold_arr)
            mf, vf = float(np.mean(fold_arr)), float(np.var(fold_arr))
            mc, vc = float(np.mean(cusp_arr)), float(np.var(cusp_arr))
            tail = float(np.mean(fold_arr > 0.5))
            rows.append((label, n, mf, vf, mc, vc, tail, dt))
            print(f'  {label}: N={n:4d}  '
                  f'Rf_mean={mf:.4f}  Rf_var={vf:.4f}  '
                  f'Rc_mean={mc:.4f}  Rc_var={vc:.4f}  '
                  f'P(Rf>0.5)={tail:.4f}  [{dt:.0f}s]', flush=True)
    return rows


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Decomposition analysis: isolate LOS vs subhalo effects on R_fold/R_cusp')
    parser.add_argument('--n-draws', type=int, default=2000, help='Number of lens draws to attempt')
    parser.add_argument('--parallel', action='store_true', help='Parallel processing')
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--quick', action='store_true', help='Quick mode: 500 draws')
    args = parser.parse_args()
    if args.quick:
        args.n_draws = 500
        args.parallel = True
    n_workers = args.workers or os.cpu_count() if args.parallel else 1

    f_sub_list = [0.0, 0.005, 0.020]
    configs = {
        'Smooth':       {'substructure': False, 'los': False, 'snr': 20},
        'Smooth+LOS':   {'substructure': False, 'los': True,  'snr': 20},
        'Smooth+Sub':   {'substructure': True,  'los': False, 'snr': 20},
        'Smooth+LOS+Sub': {'substructure': True,  'los': True,  'snr': 20},
    }

    print(f'Preparing {args.n_draws} lens draws (finding quads)...', flush=True)
    t0 = time.perf_counter()
    quads = prepare_macro_quads(args.n_draws, 42)
    print(f'Prep done: {len(quads)} quads in {time.perf_counter()-t0:.0f}s', flush=True)

    print(f'\n=== Decomposition Table ({len(quads)} quads, {len(f_sub_list)} f_sub values) ===')
    rows = run_decomposition_table(quads, f_sub_list, configs, n_workers)

    print('\n' + '=' * 120)
    print(f'{"Config":20s} {"N":>5s} {"Rf_mean":>8s} {"Rf_var":>8s} {"Rc_mean":>8s} {"Rc_var":>8s} {"P(Rf>0.5)":>10s} {"Time":>6s}')
    print('-' * 120)
    for label, n, mf, vf, mc, vc, tail, dt in rows:
        print(f'{label:20s} {n:5d} {mf:8.4f} {vf:8.4f} {mc:8.4f} {vc:8.4f} {tail:10.4f} {dt:6.0f}s')

    print('\n=== Analysis ===')
    smooth_row = next(r for r in rows if r[0].startswith('Smooth  '))
    los_row = next(r for r in rows if r[0].startswith('Smooth+LOS  '))
    sub_row = next(r for r in rows if r[0].startswith('Smooth+Sub  '))
    full_row = next(r for r in rows if r[0].startswith('Smooth+LOS+Sub  '))

    print(f'  Base smooth R_fold variance:    {smooth_row[3]:.6f}')
    print(f'  +LOS added variance:            {los_row[3] - smooth_row[3]:.6f} (total: {los_row[3]:.6f})')
    print(f'  +Sub (f_sub=0) added variance:   {sub_row[3] - smooth_row[3]:.6f} (total: {sub_row[3]:.6f})')
    print(f'  +LOS+Sub (f_sub=0) variance:     {full_row[3]:.6f}')

    for fs, label in zip(f_sub_list, ['Smooth+Sub', 'Smooth+LOS+Sub']):
        r = next(ro for ro in rows if ro[0].startswith(f'{label}  f_sub={fs}'))
        r0 = next(ro for ro in rows if ro[0].startswith(f'{label}  f_sub=0.0'))
        print(f'  {label} f_sub={fs}: Rf_mean={r[2]:.4f} (delta={r[2]-r0[2]:.4f})  tail={r[6]:.4f}')


if __name__ == '__main__':
    main()
