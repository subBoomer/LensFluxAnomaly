import argparse, sys, time, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))

from src.catalog_utils import build_unified_catalog
from src.compute_rmin import compute_rmin
from src.perturbation_model import SimplePerturbationModel
from scipy.stats import ks_2samp, anderson_ksamp


def tail_ratio(obs, sim, threshold=0.2):
    return float(np.mean(obs > threshold) / max(np.mean(sim > threshold), 1e-10))


def main():
    parser = argparse.ArgumentParser(description='R_min analysis: compare flux-ratio anomalies against LCDM')
    parser.add_argument('--n-sim', type=int, default=50000, help='Number of simulated systems')
    parser.add_argument('--delta-r', type=float, default=0.2, help='Radial pairing threshold (arcsec)')
    parser.add_argument('--mode', choices=['radio', 'optical', 'all'], default='all',
                        help='Which sample to use')
    parser.add_argument('--epsilon', type=float, default=3e-9, help='Perturbation scaling')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    inc_radio = args.mode in ('radio', 'all')
    inc_optical = args.mode in ('optical', 'all')
    systems = build_unified_catalog(include_radio=inc_radio, include_optical=inc_optical,
                                    deduplicate=(args.mode == 'all'))
    print(f'Catalog: {len(systems)} systems (mode={args.mode})')

    obs_results = []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes'], args.delta_r)
        if rmin is not None:
            obs_results.append(rmin)
            print(f'  {s["name"]:15s}  band={s["band"]:20s}  R_min = {rmin:.4f}')
        else:
            print(f'  {s["name"]:15s}  R_min = None')
    obs_arr = np.array(obs_results)
    print(f'  ---> {len(obs_arr)}/{len(systems)} valid R_min values')

    if len(obs_arr) < 3:
        print('Too few observed systems.')
        sys.exit(1)

    print(f'\nSimulating {args.n_sim} LCDM realisations (eps={args.epsilon:.0e})...')
    t0 = time.perf_counter()
    model = SimplePerturbationModel(epsilon=args.epsilon)
    sim_arr = model.sample(args.n_sim, seed=args.seed)
    dt = time.perf_counter() - t0

    ks = ks_2samp(obs_arr, sim_arr)
    try:
        ad = anderson_ksamp([obs_arr.tolist(), sim_arr.tolist()])
        ad_sig = ad.significance_level
    except Exception:
        ad_sig = None
    tail_obs = np.mean(obs_arr > 0.2)
    tail_sim = np.mean(sim_arr > 0.2)
    T = tail_ratio(obs_arr, sim_arr)

    print('\n' + '=' * 60)
    print(f'  Mode:                   {args.mode}')
    print(f'  Systems:                {len(obs_arr)}')
    print(f'  Simulations:            {len(sim_arr)}')
    print(f'  Delta r:                {args.delta_r} arcsec')
    print(f'')
    print(f'  Observed R_min: mean={np.mean(obs_arr):.4f}  std={np.std(obs_arr):.4f}')
    print(f'    values: {", ".join(f"{r:.4f}" for r in sorted(obs_arr))}')
    print(f'  Simulated R_min: mean={np.mean(sim_arr):.4f}  std={np.std(sim_arr):.4f}')
    print(f'')
    print(f'  KS test:                D={ks.statistic:.4f}  p={ks.pvalue:.4f}')
    if ad_sig is not None:
        print(f'  Anderson-Darling:        sig={ad_sig:.4f}')
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
            ax1.axvline(r, color='red', alpha=0.7, linewidth=1.2)
        ax1.set_xlabel('R_min')
        ax1.set_ylabel('Density')
        ax1.set_title(f'R_min distribution ({args.mode}, N={len(obs_arr)})')
        ax1.legend()
        obs_sorted = np.sort(obs_arr)
        obs_ecdf_y = np.arange(1, len(obs_sorted) + 1) / len(obs_sorted)
        sim_sorted = np.sort(sim_arr)
        sim_ecdf_y = np.arange(1, len(sim_sorted) + 1) / len(sim_sorted)
        ax2.step(obs_sorted, obs_ecdf_y, where='post', color='red', linewidth=2, label=f'Obs (N={len(obs_arr)})')
        ax2.step(sim_sorted, sim_ecdf_y, where='post', color='gray', linewidth=1.5, label=f'LCDM sim')
        ax2.set_xlabel('R_min')
        ax2.set_ylabel('CDF')
        ax2.set_title(f'AD sig={ad_sig:.3f}, T={T:.2f}' if ad_sig else f'KS p={ks.pvalue:.3f}, T={T:.2f}')
        ax2.legend()
        plt.tight_layout()
        out = Path(f'outputs/rmin_{args.mode}.png')
        fig.savefig(str(out), dpi=150)
        print(f'\n  Plot: {out}')
        plt.close(fig)
    except Exception as e:
        print(f'\n  Plot skipped: {e}')


if __name__ == '__main__':
    main()
