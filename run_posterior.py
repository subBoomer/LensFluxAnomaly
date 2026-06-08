import numpy as np
from pathlib import Path
from data.radio_quads import build_observed_catalog as get_obs
from src.statistic import compute_rfold
from src.inference import run_inference

outputs_dir = Path('outputs')

R_obs = []
for s in get_obs():
    r = compute_rfold(np.array(s['x_arcsec']), np.array(s['y_arcsec']),
        np.array(s['parity']), np.array(s['fluxes_mjy']))
    if r >= 0:
        R_obs.append(r)
R_obs = np.array(R_obs)
print(f'Observed: {len(R_obs)}')

f_grid = []
R_grid = []
for fs in [0.0, 0.001, 0.002, 0.005, 0.01, 0.02]:
    ckpt = outputs_dir / f'fsub_{fs:.3f}.npy'
    if ckpt.exists():
        R = np.load(ckpt)
        f_grid.append(fs)
        R_grid.append(R)
        print(f'  f_sub={fs:.4f}: N={len(R)}, mean={np.mean(R):.4f}')

print()
post = run_inference(f_grid, R_grid, R_obs, prior='uniform')
print(f'Best-fit f_sub = {post["best_fit"]:.4f}')
print(f'Median   f_sub = {post["median"]:.4f}')
print(f'68% CI        = [{post["ci_68_low"]:.4f}, {post["ci_68_high"]:.4f}]')
print(f'Mean +- std   = {post["mean"]:.4f} +- {post["std"]:.4f}')

np.savez(outputs_dir / 'fsub_posterior.npz',
    f_sub_vals=post['f_sub_vals'], posterior=post['posterior'],
    best_fit=np.array([post['best_fit']]),
    median=np.array([post['median']]),
    ci_68_low=np.array([post['ci_68_low']]),
    ci_68_high=np.array([post['ci_68_high']]))

# quick plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.fill_between(post['f_sub_vals'], post['posterior'], alpha=0.3)
    ax1.plot(post['f_sub_vals'], post['posterior'])
    ax1.axvline(post['best_fit'], color='C1', ls='--', label=f'best={post["best_fit"]:.4f}')
    ax1.axvline(post['ci_68_low'], color='gray', ls=':', label=f'68% CI')
    ax1.axvline(post['ci_68_high'], color='gray', ls=':')
    ax1.set(xlabel='f_sub', ylabel='posterior density', title='f_sub posterior')
    ax1.legend(fontsize=8)
    for fs, R in zip(f_grid, R_grid):
        label = f'f_sub={fs:.3f}'
        ax2.hist(R, bins=15, histtype='step', density=True, label=label, alpha=0.7)
    ax2.set(xlabel='R_fold', ylabel='density', title='Simulated R_fold distributions')
    ax2.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(str(outputs_dir / 'fsub_posterior.png'), dpi=150)
    print('Saved fsub_posterior.png')
except Exception as e:
    print(f'Plot skip: {e}')
