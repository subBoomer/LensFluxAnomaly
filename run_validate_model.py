import argparse, os, sys, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))

from src.catalog_utils import build_unified_catalog
from src.compute_rmin import compute_rmin
from src.perturbation_model import SimplePerturbationModel
from src.population import LensPopulation
from src.lens_model import MacroLens
from src.substructure import SubhaloPopulation, LOSPopulation
from src.selection import passes_selection
from src.noise_model import RadioNoise
from scipy.stats import ks_2samp, anderson_ksamp


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


def simulate_one_lenstronomy(seed, f_sub):
    rng = np.random.default_rng(seed)
    pop = LensPopulation(rng)
    sub_pop = SubhaloPopulation(f_sub=f_sub)
    los_pop = LOSPopulation()
    radio_noise = RadioNoise()

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
        return None
    if not passes_selection(
        macro_result['theta_x'], macro_result['theta_y'],
        macro_result['mu'], theta['source_flux'],
    ):
        return None
    subhalos = sub_pop.realise(theta_E, theta['z_l'], rng)
    macro.add_substructure(subhalos)
    full_result = macro.solve(bx, by, num_random=2, search_window=2)
    if full_result is None or full_result['n_images'] < 4:
        return None
    F_true = np.abs(full_result['mu']) * theta['source_flux']
    F_obs = radio_noise.apply(F_true, rng)
    x = full_result['theta_x']
    y = full_result['theta_y']
    return x, y, F_obs


def simulate_batch_lenstronomy(seeds, f_sub):
    results = []
    for seed in seeds:
        out = simulate_one_lenstronomy(seed, f_sub)
        if out is not None:
            results.append(out)
    return results


def tail_ratio(obs, sim, threshold=0.2):
    return float(np.mean(obs > threshold) / max(np.mean(sim > threshold), 1e-10))


def run_comparison(obs_arr, sim_arr, label):
    if len(obs_arr) < 3 or len(sim_arr) < 10:
        print(f'  {label}: too few samples ({len(obs_arr)} obs, {len(sim_arr)} sim)')
        return
    ks = ks_2samp(obs_arr, sim_arr)
    try:
        ad = anderson_ksamp([obs_arr.tolist(), sim_arr.tolist()])
        ad_sig = ad.significance_level
    except Exception:
        ad_sig = None
    T = tail_ratio(obs_arr, sim_arr)
    print(f'  {label}:')
    print(f'    Obs: mean={np.mean(obs_arr):.4f} std={np.std(obs_arr):.4f} (N={len(obs_arr)})')
    print(f'    Sim: mean={np.mean(sim_arr):.4f} std={np.std(sim_arr):.4f} (N={len(sim_arr)})')
    print(f'    KS: D={ks.statistic:.4f} p={ks.pvalue:.4f}')
    if ad_sig is not None:
        print(f'    AD: sig={ad_sig:.4f}')
    print(f'    Tail ratio T: {T:.3f}')
    return {'ks': ks, 'ad_sig': ad_sig, 'T': T, 'obs': obs_arr, 'sim': sim_arr}


def main():
    parser = argparse.ArgumentParser(description='Phase B: cross-validate R_min with full lenstronomy pipeline')
    parser.add_argument('--n-sim', type=int, default=2000, help='Number of lenstronomy realisations')
    parser.add_argument('--n-perturb', type=int, default=50000, help='Number of simple perturbation realisations')
    parser.add_argument('--f-sub', type=float, default=0.005, help='Subhalo mass fraction (DMO=0.005)')
    parser.add_argument('--delta-r', type=float, default=0.2, help='Radial pairing threshold (arcsec)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--parallel', action='store_true', help='Use multiprocessing')
    parser.add_argument('--quick', action='store_true', help='Quick test: 100 draws')
    parser.add_argument('--force', action='store_true', help='Re-run even if results exist')
    args = parser.parse_args()

    if args.quick:
        args.n_sim = 100
        args.n_perturb = 10000

    print('=' * 60)
    print('Phase B: lenstronomy cross-validation for R_min')
    print('=' * 60)

    systems = build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True)
    obs_results = []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes'], args.delta_r)
        if rmin is not None:
            obs_results.append(rmin)
    obs_arr = np.array(obs_results)
    print(f'\nObserved catalog: {len(systems)} systems, {len(obs_arr)} valid R_min')

    print(f'\nGenerating {args.n_sim} lenstronomy realisations (f_sub={args.f_sub})...')
    t0 = time.perf_counter()
    seeds = list(range(args.seed, args.seed + args.n_sim))
    lens_results = []
    if args.parallel:
        n_workers = os.cpu_count()
        batch_size = max(1, args.n_sim // (n_workers * 4))
        batches = [seeds[i:i + batch_size] for i in range(0, len(seeds), batch_size)]
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            fut_to_idx = {pool.submit(simulate_batch_lenstronomy, b, args.f_sub): i
                          for i, b in enumerate(batches)}
            done = 0
            for fut in as_completed(fut_to_idx):
                batch_results = fut.result()
                lens_results.extend(batch_results)
                done += 1
                if done % 10 == 0 or done == len(batches):
                    print(f'  Batch {done}/{len(batches)} ({len(lens_results)} quads so far)', flush=True)
    else:
        for i, seed in enumerate(seeds):
            out = simulate_one_lenstronomy(seed, args.f_sub)
            if out is not None:
                lens_results.append(out)
            if (i + 1) % 100 == 0:
                print(f'  {i+1}/{args.n_sim} ({len(lens_results)} quads)', flush=True)
    dt = time.perf_counter() - t0

    lens_rmin = []
    for x, y, F in lens_results:
        r = compute_rmin(x, y, F, args.delta_r)
        if r is not None:
            lens_rmin.append(r)
    lens_arr = np.array(lens_rmin)
    print(f'  Done in {dt:.1f}s: {len(lens_results)} quads, {len(lens_arr)} valid R_min')

    print(f'\nGenerating {args.n_perturb} simple perturbation model realisations...')
    t0 = time.perf_counter()
    model = SimplePerturbationModel()
    sim_arr = model.sample(args.n_perturb, seed=args.seed + 9999)
    print(f'  Done in {time.perf_counter() - t0:.1f}s: {len(sim_arr)} valid')

    print('\n' + '=' * 60)
    print('Comparisons')
    print('=' * 60)

    c1 = run_comparison(obs_arr, lens_arr, '1. Observed vs lenstronomy (SIE+TNFW+LOS)')
    c2 = run_comparison(obs_arr, sim_arr, '2. Observed vs simple perturbation')
    c3 = run_comparison(lens_arr, sim_arr, '3. Lenstronomy vs simple perturbation')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        bins = np.linspace(0, 0.6, 61)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        for ax, sim, obs, title, lbl in [
            (axes[0], lens_arr, obs_arr, f'Lenstronomy vs Obs\nN_obs={len(obs_arr)} N_sim={len(lens_arr)}', 'Observed'),
            (axes[1], sim_arr, obs_arr, f'Simple vs Obs\nN_obs={len(obs_arr)} N_sim={len(sim_arr)}', 'Observed'),
            (axes[2], sim_arr, lens_arr, f'Lenstronomy vs Simple\nN_lens={len(lens_arr)} N_sim={len(sim_arr)}', 'Lenstronomy'),
        ]:
            ax.hist(sim, bins=bins, density=True, alpha=0.6, color='gray', label='Simulated')
            for r in obs:
                ax.axvline(r, color='red', alpha=0.7, linewidth=1.2, label=lbl if r == obs[0] else '')
            ax.set_xlabel('R_min')
            ax.set_ylabel('Density')
            ax.set_title(title)
            ax.legend()

        plt.tight_layout()
        out = Path('outputs/phase_b_validation.png')
        fig.savefig(str(out), dpi=150)
        print(f'\nPlot saved: {out}')
        plt.close(fig)
    except Exception as e:
        print(f'\nPlot skipped: {e}')

    result_path = Path('outputs/phase_b_validation.npz')
    np.savez(result_path,
             obs_rmin=obs_arr,
             lenstronomy_rmin=lens_arr,
             simple_rmin=sim_arr,
             )
    print(f'Results saved: {result_path}')


if __name__ == '__main__':
    main()
