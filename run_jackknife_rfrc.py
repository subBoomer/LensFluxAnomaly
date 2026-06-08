import sys, warnings, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
warnings.filterwarnings('ignore', category=UserWarning)

from src.catalog_utils import build_unified_catalog
from src.statistic import compute_rfold, compute_cusp_statistic
from src.perturbation_model import SimplePerturbationModel
from scipy.stats import anderson_ksamp

N_SIM = 50000
N_PERM = 500
SEED = 42


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


def main():
    print('=' * 60)
    print('Jackknife: R_fold and R_cusp')
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
            obs_entries.append({'name': s['name'], 'rfold': rfold, 'rcusp': rcusp})
    obs_rfold = np.array([e['rfold'] for e in obs_entries])
    obs_rcusp = np.array([e['rcusp'] for e in obs_entries])
    print(f'\nObserved: N={len(obs_entries)} radio systems with parity')

    data = np.load('outputs/rfold_rcusp_analysis.npz', allow_pickle=True)
    null_rfold = data['cdm_rfold']
    null_rcusp = data['cdm_rcusp']
    print(f'\nNull: R_fold N={len(null_rfold)}, R_cusp N={len(null_rcusp)}')
    print(f'  R_fold: mean={null_rfold.mean():.4f} std={null_rfold.std():.4f}')
    print(f'  R_cusp: mean={null_rcusp.mean():.4f} std={null_rcusp.std():.4f}')

    for stat, obs_arr, null_arr, label in [
        ('rfold', obs_rfold, null_rfold, 'R_fold'),
        ('rcusp', obs_rcusp, null_rcusp, 'R_cusp'),
    ]:
        ad_full = anderson_ksamp([obs_arr.tolist(), null_arr.tolist()])
        p_full = permutation_ad_pvalue(obs_arr, null_arr)
        print(f'\n{"="*60}')
        print(f'  {label} — Full sample')
        print(f'{"="*60}')
        print(f'  Full: N={len(obs_arr)} AD_stat={ad_full.statistic:.2f} sig={ad_full.significance_level:.4f} p_perm={p_full:.4f}')
        print()
        print(f'  {"System":18s} {"Value":8s} {"AD_stat":8s} {"p_perm":8s} {"Delta":8s}')
        print(f'  {"-"*50}')
        for i, entry in enumerate(obs_entries):
            reduced = np.delete(obs_arr, i)
            ad = anderson_ksamp([reduced.tolist(), null_arr.tolist()])
            p = permutation_ad_pvalue(reduced, null_arr)
            delta = ad.statistic - ad_full.statistic
            p_str = f'{p:.4f}' if p >= 0.0001 else '<.0001'
            val = entry[stat]
            print(f'  {entry["name"]:18s} {val:.4f}   {ad.statistic:8.2f} {p_str:8s} {delta:+7.2f}')
        print(f'  {"-"*50}')
        print(f'  AD stat range: {min(ad_full.statistic + delta for delta in [0]):.2f} (reference)')


if __name__ == '__main__':
    main()
