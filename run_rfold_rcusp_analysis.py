import sys, time, warnings, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
warnings.filterwarnings('ignore', category=UserWarning)

from src.catalog_utils import build_unified_catalog
from src.statistic import compute_rfold, compute_cusp_statistic
from src.perturbation_model import SimplePerturbationModel
from scipy.stats import anderson_ksamp, ks_2samp

N_SIM = 50000
N_PERM = 2000
SEED = 42


def tail_frac(arr, threshold=0.2):
    return np.mean(arr > threshold) if len(arr) > 0 else 0.0


def permutation_ad_pvalue(obs, null, n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    combined = np.concatenate([obs, null])
    n_obs = len(obs)
    ad_obs = anderson_ksamp([obs.tolist(), null.tolist()]).statistic
    count = 1
    for _ in range(n_perm):
        rng.shuffle(combined)
        try:
            ad_p = anderson_ksamp([combined[:n_obs].tolist(), combined[n_obs:].tolist()]).statistic
        except Exception:
            ad_p = 0
        if ad_p >= ad_obs:
            count += 1
    return count / (n_perm + 1)


def compute_rfold_rcusp_from_abstract(F_obs, x_img, y_img):
    parity = np.array([1, -1, 1, -1])
    rfold = compute_rfold(x_img, y_img, parity.astype(float), F_obs)
    rcusp = compute_cusp_statistic(x_img, y_img, parity, F_obs)
    return rfold, rcusp


def sample_perturbation_stats(model, n_sim, seed=SEED, verbose=True):
    rng = np.random.default_rng(seed)
    rmin_list, rfold_list, rcusp_list = [], [], []
    for i in range(n_sim):
        result = model.realise(seed=rng.integers(0, 2**31), return_all=True)
        if result is None or result[0] is None:
            continue
        rmin, x, y, F, *_ = result
        rmin_list.append(rmin)
        rf, rc = compute_rfold_rcusp_from_abstract(F, x, y)
        if rf >= 0:
            rfold_list.append(rf)
        if rc >= 0:
            rcusp_list.append(rc)
        if verbose and (i + 1) % 10000 == 0:
            print(f'  {i+1}/{n_sim} ({len(rmin_list)} valid)', flush=True)
    return np.array(rmin_list), np.array(rfold_list), np.array(rcusp_list)


def run_comparison(label, obs, sim, null_tail_frac=None):
    if len(obs) < 3 or len(sim) < 10:
        print(f'  {label}: too few samples')
        return None
    ad = anderson_ksamp([obs.tolist(), sim.tolist()])
    p = permutation_ad_pvalue(obs, sim)
    ks = ks_2samp(obs, sim)
    T = tail_frac(obs) / max(tail_frac(sim), 1e-10)
    print(f'  {label}:')
    print(f'    Obs: mean={np.mean(obs):.4f} std={np.std(obs):.4f} (N={len(obs)})')
    print(f'    Sim: mean={np.mean(sim):.4f} std={np.std(sim):.4f} (N={len(sim)})')
    print(f'    KS:  D={ks.statistic:.4f} p={ks.pvalue:.4f}')
    print(f'    AD:  stat={ad.statistic:.2f} sig={ad.significance_level:.4f} p_perm={p:.4f}')
    print(f'    Tail T={T:.2f}')
    return {'ad_stat': ad.statistic, 'ad_sig': ad.significance_level, 'p_perm': p,
            'ks_stat': ks.statistic, 'ks_p': ks.pvalue, 'T': T,
            'obs_mean': np.mean(obs), 'obs_std': np.std(obs), 'n_obs': len(obs),
            'sim_mean': np.mean(sim), 'sim_std': np.std(sim), 'n_sim': len(sim)}


def main():
    print('=' * 60)
    print('R_fold and R_cusp: Full Treatment')
    print('=' * 60)

    systems = build_unified_catalog(include_radio=True, include_optical=False, deduplicate=False)
    obs_entries = []
    for s in systems:
        if s.get('parity') is None:
            continue
        x = np.array(s['x_arcsec'], dtype=float)
        y = np.array(s['y_arcsec'], dtype=float)
        F = np.array(s['fluxes'], dtype=float)
        p = np.array(s['parity'], dtype=float)
        rfold = compute_rfold(x, y, p, F)
        rcusp = compute_cusp_statistic(x, y, p, F)
        if rfold >= 0 and rcusp >= 0:
            obs_entries.append({'name': s['name'], 'rfold': rfold, 'rcusp': rcusp, 'band': s['band']})
    obs_rfold = np.array([e['rfold'] for e in obs_entries])
    obs_rcusp = np.array([e['rcusp'] for e in obs_entries])
    print(f'\nObserved radio systems with parity: {len(obs_entries)}')
    print(f'  {"System":18s} {"Band":16s} {"R_fold":8s} {"R_cusp":8s}')
    print('  ' + '-' * 50)
    for e in obs_entries:
        print(f'  {e["name"]:18s} {e["band"]:16s} {e["rfold"]:.4f}   {e["rcusp"]:.4f}')

    print(f'\nGenerating {N_SIM} CDM simulations...')
    t0 = time.perf_counter()
    model = SimplePerturbationModel()
    cdm_rmin, cdm_rfold, cdm_rcusp = sample_perturbation_stats(model, N_SIM)
    print(f'  Done in {time.perf_counter()-t0:.1f}s')
    print(f'  CDM: R_fold N={len(cdm_rfold)}, R_cusp N={len(cdm_rcusp)}')

    results = {}
    print('\n' + '=' * 60)
    print('R_fold comparisons')
    print('=' * 60)
    results['rfold'] = run_comparison('R_fold (obs vs CDM)', obs_rfold, cdm_rfold)

    print('\n' + '=' * 60)
    print('R_cusp comparisons')
    print('=' * 60)
    results['rcusp'] = run_comparison('R_cusp (obs vs CDM)', obs_rcusp, cdm_rcusp)

    results['obs_rfold'] = obs_rfold
    results['obs_rcusp'] = obs_rcusp
    results['cdm_rfold'] = cdm_rfold
    results['cdm_rcusp'] = cdm_rcusp
    results['obs_entries'] = obs_entries

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        bins = np.linspace(0, 1.0, 41)
        for ax, obs, sim, stat, label in [
            (axes[0], obs_rfold, cdm_rfold, 'R_fold', 'R_fold'),
            (axes[1], obs_rcusp, cdm_rcusp, 'R_cusp', 'R_cusp'),
        ]:
            ax.hist(sim, bins=bins, density=True, alpha=0.6, color='gray', label='CDM sim')
            for r in obs:
                ax.axvline(r, color='red', alpha=0.7, linewidth=1.2)
            ax.set_xlabel(label)
            ax.set_ylabel('Density')
            ax.set_title(f'{label}\nN_obs={len(obs)} N_sim={len(sim)}')
            ax.legend()
        plt.tight_layout()
        out = Path('outputs/rfold_rcusp_analysis.png')
        fig.savefig(str(out), dpi=150)
        print(f'\nPlot: {out}')
        plt.close(fig)
    except Exception as e:
        print(f'\nPlot skipped: {e}')

    out_path = Path('outputs/rfold_rcusp_analysis.npz')
    np.savez(out_path,
             obs_rfold=obs_rfold, obs_rcusp=obs_rcusp,
             cdm_rfold=cdm_rfold, cdm_rcusp=cdm_rcusp,
             systems=[e['name'] for e in obs_entries],
             rfold_per_sys=np.array([e['rfold'] for e in obs_entries]),
             rcusp_per_sys=np.array([e['rcusp'] for e in obs_entries]),
             )
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    main()
