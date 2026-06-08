import sys, warnings, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))

from src.catalog_utils import build_unified_catalog
from src.compute_rmin import compute_rmin
from scipy.stats import anderson_ksamp, ks_2samp

warnings.filterwarnings('ignore', category=UserWarning)


def permutation_ad_pvalue(obs, null, n_perm=5000, seed=42):
    rng = np.random.default_rng(seed)
    combined = np.concatenate([obs, null])
    n_obs = len(obs)
    ad_obs = anderson_ksamp([obs.tolist(), null.tolist()]).statistic
    count = 1
    for _ in range(n_perm):
        rng.shuffle(combined)
        obs_p = combined[:n_obs]
        null_p = combined[n_obs:]
        try:
            ad_p = anderson_ksamp([obs_p.tolist(), null_p.tolist()]).statistic
        except Exception:
            ad_p = 0
        if ad_p >= ad_obs:
            count += 1
    return count / (n_perm + 1)


def tail_frac(arr, threshold=0.2):
    return np.mean(arr > threshold)


def main():
    delta_r = 0.2

    data = np.load('outputs/phase_b_validation.npz')
    null_rmin = data['lenstronomy_rmin']
    null_tail = tail_frac(null_rmin)

    systems = build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True)
    obs_entries = []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes'], delta_r)
        if rmin is not None:
            stype = 'optical' if s['survey'] == 'LITERATURE' else 'radio'
            obs_entries.append({'name': s['name'], 'type': stype, 'rmin': rmin})
    obs_arr = np.array([e['rmin'] for e in obs_entries])

    ad_full = anderson_ksamp([obs_arr.tolist(), null_rmin.tolist()])
    p_full = permutation_ad_pvalue(obs_arr, null_rmin)
    T_full = tail_frac(obs_arr) / null_tail

    print(f'Null: mean={null_rmin.mean():.4f} tail={null_tail:.4f} (N={len(null_rmin)})')
    print(f'Full: N={len(obs_arr)} AD_stat={ad_full.statistic:.2f} AD_sig={ad_full.significance_level:.4f} p_perm={p_full:.4f} T={T_full:.1f}')
    radio_only = np.array([e['rmin'] for e in obs_entries if e['type'] == 'radio'])
    optical_only = np.array([e['rmin'] for e in obs_entries if e['type'] == 'optical'])
    print(f'Radio: N={len(radio_only)} mean={radio_only.mean():.4f} T={tail_frac(radio_only)/null_tail:.1f}')
    print(f'Opt:   N={len(optical_only)} mean={optical_only.mean():.4f} T={tail_frac(optical_only)/null_tail:.1f}')
    print()
    print(f'{"System":18s} {"Typ":4s} {"R_min":7s} {"AD_stat":8s} {"p_perm":7s} {"Tail_T":6s} {"Removed":8s}')
    print('-' * 65)

    results = []
    for i, entry in enumerate(obs_entries):
        reduced = np.delete(obs_arr, i)
        ad = anderson_ksamp([reduced.tolist(), null_rmin.tolist()])
        p = permutation_ad_pvalue(reduced, null_rmin)
        T = tail_frac(reduced) / null_tail
        p_str = f'{p:.4f}' if p >= 0.0001 else '<.0001'
        removed = 'TAIL' if entry['rmin'] > 0.2 else 'tail' if entry['rmin'] > 0.15 else 'body'
        print(f'{entry["name"]:18s} {entry["type"]:4s} {entry["rmin"]:.4f}  {ad.statistic:8.2f} {p_str:7s} {T:6.1f} {removed:8s}')
        results.append({'name': entry['name'], 'type': entry['type'], 'rmin': entry['rmin'],
                        'ad_stat': ad.statistic, 'p': p, 'T': T})

    print()
    print(f'AD stat range across jackknife: {min(r["ad_stat"] for r in results):.2f} - {max(r["ad_stat"] for r in results):.2f}')
    print(f'Reference (full): {ad_full.statistic:.2f}')
    print(f'Drivers: {", ".join(r["name"] for r in results if abs(r["ad_stat"] - ad_full.statistic) > 1.5)}')
    print(f'Survivors: {", ".join(r["name"] for r in results if abs(r["ad_stat"] - ad_full.statistic) <= 1.5)}')


if __name__ == '__main__':
    main()
