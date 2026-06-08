import sys, time, warnings, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
warnings.filterwarnings('ignore', category=UserWarning)

from src.catalog_utils import build_unified_catalog
from src.compute_rmin import compute_rmin
from src.statistic import compute_rfold, compute_cusp_statistic
from src.perturbation_model import SimplePerturbationModel
from src.wdm_model import WDMPerturbationModel
from scipy.stats import anderson_ksamp

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
    """Infer parity from abstract quad geometry and compute R_fold, R_cusp."""
    parity = np.array([1, -1, 1, -1])
    rfold = compute_rfold(x_img, y_img, parity.astype(float), F_obs)
    rcusp = compute_cusp_statistic(x_img, y_img, parity, F_obs)
    return rfold, rcusp


def sample_perturbation_stats(model, n_sim, seed=42, verbose=True):
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


def run_comparison(label, obs, sim, null_tail):
    if len(obs) < 3 or len(sim) < 10:
        print(f'  {label}: too few samples')
        return
    ad = anderson_ksamp([obs.tolist(), sim.tolist()])
    p = permutation_ad_pvalue(obs, sim)
    T = tail_frac(obs) / max(tail_frac(sim), 1e-10)
    print(f'  {label}:')
    print(f'    Obs: mean={np.mean(obs):.4f} std={np.std(obs):.4f}  Sim: mean={np.mean(sim):.4f} std={np.std(sim):.4f}')
    print(f'    AD stat={ad.statistic:.2f}  AD sig={ad.significance_level:.4f}  p_perm={p:.4f}  T={T:.2f}')


def main():
    print('=' * 60)
    print('Phase F: WDM predictions + R_fold/R_cusp full treatment')
    print('=' * 60)

    systems = build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True)
    obs_rmin, obs_rfold, obs_rcusp = [], [], []
    for s in systems:
        rmin = compute_rmin(s['x_arcsec'], s['y_arcsec'], s['fluxes'], 0.2)
        if rmin is not None:
            obs_rmin.append(rmin)
        if 'parity' in s and s['parity'] is not None and len(s['parity']) == 4:
            rfold = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                                  np.array(s['parity'], dtype=float), np.array(s['fluxes']))
            rcusp = compute_cusp_statistic(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
                                           np.array(s['parity']), np.array(s['fluxes']))
            if rfold >= 0:
                obs_rfold.append(rfold)
            if rcusp >= 0:
                obs_rcusp.append(rcusp)
    obs_rmin = np.array(obs_rmin)
    obs_rfold = np.array(obs_rfold)
    obs_rcusp = np.array(obs_rcusp)
    print(f'\nObserved: R_min N={len(obs_rmin)}, R_fold N={len(obs_rfold)}, R_cusp N={len(obs_rcusp)}')

    print(f'\nGenerating {N_SIM} CDM realisations...')
    t0 = time.perf_counter()
    model_cdm = SimplePerturbationModel()
    cdm_rmin, cdm_rfold, cdm_rcusp = sample_perturbation_stats(model_cdm, N_SIM)
    print(f'  Done in {time.perf_counter()-t0:.1f}s')
    print(f'  CDM: R_min N={len(cdm_rmin)}, R_fold N={len(cdm_rfold)}, R_cusp N={len(cdm_rcusp)}')

    print('\n' + '-' * 60)
    print('Part 1: WDM predictions')
    print('-' * 60)
    for m_wdm in [3, 5, 7]:
        print(f'\n  WDM m_x={m_wdm} keV (m_hm={WDMPerturbationModel(m_wdm_keV=m_wdm).m_hm:.1e} Msun)')
        t0 = time.perf_counter()
        model_wdm = WDMPerturbationModel(m_wdm_keV=m_wdm)
        wdm_rmin, wdm_rfold, wdm_rcusp = sample_perturbation_stats(model_wdm, N_SIM)
        print(f'  Done in {time.perf_counter()-t0:.1f}s')
        print(f'  WDM: R_min N={len(wdm_rmin)}, R_fold N={len(wdm_rfold)}, R_cusp N={len(wdm_rcusp)}')
        run_comparison('WDM vs CDM (R_min)', wdm_rmin, cdm_rmin, tail_frac(cdm_rmin))
        print()

    print('-' * 60)
    print('Part 2: CDM R_fold and R_cusp vs observed')
    print('-' * 60)
    run_comparison('R_min (CDM vs obs)', obs_rmin, cdm_rmin, tail_frac(cdm_rmin))
    run_comparison('R_fold (CDM vs obs)', obs_rfold, cdm_rfold, tail_frac(cdm_rfold))
    run_comparison('R_cusp (CDM vs obs)', obs_rcusp, cdm_rcusp, tail_frac(cdm_rcusp))


if __name__ == '__main__':
    main()
