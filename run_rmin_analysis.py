import argparse, sys, time, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))

from data.radio_quads import build_observed_catalog, RADIO_QUAD_CATALOG
from src.compute_rmin import compute_rmin
from src.perturbation_model import SimplePerturbationModel
from scipy.stats import ks_2samp, anderson_ksamp


def tail_ratio(obs, sim, threshold=0.2):
    T = np.mean(obs > threshold) / max(np.mean(sim > threshold), 1e-10)
    return float(T)


def main():
    parser = argparse.ArgumentParser(description='R_min analysis: compare flux-ratio anomalies against LCDM')
    parser.add_argument('--n-sim', type=int, default=50000, help='Number of simulated systems')
    parser.add_argument('--delta-r', type=float, default=0.2, help='Radial pairing threshold (arcsec)')
    parser.add_argument('--include-optical', action='store_true', help='Include optical systems (PG1115)')
    parser.add_argument('--epsilon', type=float, default=5e-9, help='Perturbation scaling')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    if args.include_optical:
        systems = [s for s in RADIO_QUAD_CATALOG if not s.get("placeholder", False)]
        systems = [s for s in systems if s['z_l'] is not None]
    else:
        systems = build_observed_catalog()
    print(f'Observed systems: {len(systems)}')

    obs_results = []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes_mjy'], args.delta_r)
        if rmin is not None:
            obs_results.append(rmin)
            print(f'  {s["name"]:15s}  R_min = {rmin:.4f}  (N_r={len(s["x_arcsec"])})')
        else:
            print(f'  {s["name"]:15s}  R_min = None  (no paired images within dr={args.delta_r})')
    obs_arr = np.array(obs_results)
    print(f'  ---> {len(obs_arr)}/{len(systems)} systems have valid R_min')

    if len(obs_arr) < 3:
        print('Too few observed systems -- cannot proceed.')
        sys.exit(1)

    print(f'\nSimulating {args.n_sim} LCDM realisations with perturbation model...')
    t0 = time.perf_counter()
    model = SimplePerturbationModel(epsilon=args.epsilon)
    sim_arr = model.sample(args.n_sim, seed=args.seed)
    dt = time.perf_counter() - t0
    print(f'  Done: {len(sim_arr)}/{args.n_sim} valid quads in {dt:.0f}s ({len(sim_arr)/dt:.0f}/s)')

    ks = ks_2samp(obs_arr, sim_arr)
    try:
        ad = anderson_ksamp([obs_arr.tolist(), sim_arr.tolist()])
        ad_stat, ad_crit, ad_sig = ad.statistic, ad.critical_values, ad.significance_level
    except Exception:
        ad_stat, ad_crit, ad_sig = None, None, None
    tail_obs = np.mean(obs_arr > 0.2)
    tail_sim = np.mean(sim_arr > 0.2)
    T = tail_ratio(obs_arr, sim_arr)

    print('\n' + '=' * 60)
    print(f'  Systems:                {len(systems)} (observed), {len(obs_arr)} valid')
    print(f'  Simulations:            {len(sim_arr)}')
    print(f'  dr:                     {args.delta_r} arcsec')
    print(f'  eps:                    {args.epsilon}')
    print(f'')
    print(f'  Observed R_min: mean={np.mean(obs_arr):.4f}, std={np.std(obs_arr):.4f}')
    print(f'    per system: {", ".join(f"{r:.4f}" for r in obs_arr)}')
    print(f'  Simulated R_min: mean={np.mean(sim_arr):.4f}, std={np.std(sim_arr):.4f}')
    print(f'')
    print(f'  KS test:                D={ks.statistic:.4f}, p={ks.pvalue:.4f}')
    if ad_stat is not None:
        print(f'  Anderson-Darling:        stat={ad_stat:.4f}, sig={ad_sig:.4f}')
    print(f'')
    print(f'  Tail ratio P(R>0.2):')
    print(f'    Observed:             {tail_obs:.4f}  ({np.sum(obs_arr > 0.2)}/{len(obs_arr)})')
    print(f'    Simulated:            {tail_sim:.4f}  ({np.sum(sim_arr > 0.2)}/{len(sim_arr)})')
    print(f'    T = P_obs / P_sim:    {T:.3f}')
    print(f'  {"=" * 60}')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        bins = np.linspace(0, 0.6, 61)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.hist(sim_arr, bins=bins, density=True, alpha=0.6, label=f'LCDM sim (N={len(sim_arr)})', color='gray')
        for r in obs_arr:
            ax1.axvline(r, color='red', alpha=0.7, linewidth=1.5)
        ax1.set_xlabel('R_min')
        ax1.set_ylabel('Density')
        ax1.set_title(f'R_min distribution (dr={args.delta_r})')
        ax1.legend()
        obs_ecdf = np.sort(obs_arr)
        obs_ecdf_y = np.arange(1, len(obs_ecdf) + 1) / len(obs_ecdf)
        sim_ecdf = np.sort(sim_arr)
        sim_ecdf_y = np.arange(1, len(sim_ecdf) + 1) / len(sim_ecdf)
        ax2.step(obs_ecdf, obs_ecdf_y, where='post', color='red', linewidth=2, label=f'Obs (N={len(obs_arr)})')
        ax2.step(sim_ecdf, sim_ecdf_y, where='post', color='gray', linewidth=1.5, label=f'LCDM sim')
        ax2.set_xlabel('R_min')
        ax2.set_ylabel('CDF')
        ax2.set_title(f'KS p={ks.pvalue:.4f}, T={T:.2f}')
        ax2.legend()
        plt.tight_layout()
        out = Path('outputs/rmin_comparison.png')
        fig.savefig(str(out), dpi=150)
        print(f'\n  Plot saved: {out}')
        plt.close(fig)
    except Exception as e:
        print(f'\n  Plot skipped: {e}')


if __name__ == '__main__':
    main()
