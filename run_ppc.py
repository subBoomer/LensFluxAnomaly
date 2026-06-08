import sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from scipy.stats import ks_2samp, anderson_ksamp
N_PPC = 10000
N_OBS = 8
SEED = 42

def permutation_ad_pvalue(obs, null, n_perm=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    combined = np.concatenate([obs, null])
    n_obs = len(obs)
    ad_obs = anderson_ksamp([obs.tolist(), null.tolist()], variant='midrank', method='asymptotic').statistic
    count = 1
    for _ in range(n_perm):
        rng.shuffle(combined)
        try:
            ad_p = anderson_ksamp([combined[:n_obs].tolist(), combined[n_obs:].tolist()], variant='midrank', method='asymptotic').statistic
        except Exception:
            ad_p = 0
        if ad_p >= ad_obs:
            count += 1
    return count / (n_perm + 1)

def main():
    print('=' * 60)
    print('Posterior Predictive Check')
    print('=' * 60)
    data = np.load('outputs/rfold_rcusp_analysis.npz', allow_pickle=True)
    cdm_rfold = data['cdm_rfold']
    obs_rfold = data['obs_rfold']
    print(f'CDM simulations:  {len(cdm_rfold)}')
    print(f'Observed (radio): {len(obs_rfold)}')
    print(f'Observed R_fold:  {np.round(obs_rfold, 3)}')
    print(f'Observed mean:    {np.mean(obs_rfold):.4f}')
    print()
    rng = np.random.default_rng(SEED)
    ppc_means = []
    ppc_tails = []
    ppc_ad_stats = []
    ppc_ks_stats = []
    for i in range(N_PPC):
        sample = rng.choice(cdm_rfold, size=N_OBS, replace=True)
        ppc_means.append(np.mean(sample))
        ppc_tails.append(np.mean(sample > 0.2))
        ad = anderson_ksamp([obs_rfold.tolist(), sample.tolist()], variant='midrank', method='asymptotic')
        ppc_ad_stats.append(ad.statistic)
        ks = ks_2samp(obs_rfold, sample)
        ppc_ks_stats.append(ks.statistic)
        if (i + 1) % 2000 == 0:
            print(f'  PPC {i + 1}/{N_PPC}', flush=True)
    ppc_means = np.array(ppc_means)
    ppc_tails = np.array(ppc_tails)
    ppc_ad_stats = np.array(ppc_ad_stats)
    ppc_ks_stats = np.array(ppc_ks_stats)
    obs_mean = np.mean(obs_rfold)
    obs_tail = np.mean(obs_rfold > 0.2)
    obs_ad = anderson_ksamp([obs_rfold.tolist(), cdm_rfold.tolist()], variant='midrank', method='asymptotic').statistic
    obs_ks = ks_2samp(obs_rfold, cdm_rfold).statistic
    p_mean = np.mean(ppc_means >= obs_mean)
    p_tail = np.mean(ppc_tails >= obs_tail)
    p_ad = np.mean(ppc_ad_stats >= obs_ad)
    p_ks = np.mean(ppc_ks_stats >= obs_ks)
    print(f'\n=== Results ===')
    print(f'\n--- Test statistic: mean(R_fold) ---')
    print(f'Observed: {obs_mean:.4f}')
    print(f'PPC distribution: mean={np.mean(ppc_means):.4f}, std={np.std(ppc_means):.4f}')
    print(f'PPC p-value (P(mean >= obs)): {p_mean:.4f}')
    print(f'Significance: {('>3sigma' if p_mean < 0.003 else '>2sigma' if p_mean < 0.05 else 'not significant')}')
    print(f'\n--- Test statistic: P(R_fold > 0.2) ---')
    print(f'Observed: {obs_tail:.4f}')
    print(f'PPC p-value: {p_tail:.4f}')
    print(f'\n--- Test statistic: AD ---')
    print(f'Observed AD stat: {obs_ad:.2f}')
    print(f'PPC p-value: {p_ad:.4f}')
    ad_p_perm = permutation_ad_pvalue(obs_rfold, cdm_rfold)
    print(f'AD permutation p-value (direct): {ad_p_perm:.4f}')
    print(f'\n--- Test statistic: KS ---')
    print(f'Observed KS stat: {obs_ks:.4f}')
    print(f'PPC p-value: {p_ks:.4f}')
    print(f'\n=== Summary ===')
    print(f'{'Method':25s} {'p-value':>10s} {'Significant?':>12s}')
    print('-' * 49)
    for label, p in [('PPC (mean R_fold)', p_mean), ('PPC (tail frac)', p_tail), ('PPC (AD stat)', p_ad), ('PPC (KS stat)', p_ks), ('AD permutation', ad_p_perm)]:
        sig = 'YES' if p < 0.05 else 'no'
        print(f'{label:25s} {p:>10.4f} {sig:>12s}')
    np.savez('outputs/ppc_results.npz', ppc_means=ppc_means, ppc_tails=ppc_tails, ppc_ad_stats=ppc_ad_stats, ppc_ks_stats=ppc_ks_stats, obs_mean=obs_mean, obs_tail=obs_tail, p_mean=p_mean, p_tail=p_tail, p_ad=p_ad, p_ks=p_ks)
    print(f'\nSaved: outputs/ppc_results.npz')
if __name__ == '__main__':
    main()