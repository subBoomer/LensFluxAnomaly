import sys, numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from data.radio_quads import build_observed_catalog
from src.statistic import compute_rfold, compute_cusp_statistic
from src.comparison import compare

outputs_dir = Path(__file__).parent.parent / 'outputs'

systems = build_observed_catalog()
R_obs, C_obs = [], []
for s in systems:
    x = np.array(s['x_arcsec']); y = np.array(s['y_arcsec'])
    p = np.array(s['parity']); f = np.array(s['fluxes_mjy'])
    r = compute_rfold(x, y, p, f)
    if r >= 0:
        R_obs.append(r)
        C_obs.append(compute_cusp_statistic(x, y, p, f))
R_obs, C_obs = np.array(R_obs), np.array(C_obs)

configs = [
    ('macro_only',      'SIE+SHEAR only',           '#66c2a5'),
    ('macro_los',       '+ LOS',                     '#fc8d62'),
    ('macro_subhalos',  '+ subhalos',                '#8da0cb'),
    ('macro_sub_los',   '+ subhalos + LOS',          '#e78ac3'),
    ('macro_sub_sel',   '+ selection',               '#a6d854'),
    ('full',            'full (LOS+sub+sel+noise)',  '#ffd92f'),
]
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for stat_name, obs_vals, xlabel, title in [
    ('R_fold', R_obs, r'$R_{\mathrm{fold}}$', 'Fold statistic'),
    ('R_cusp', C_obs, r'$R_{\mathrm{cusp}}$', 'Cusp statistic'),
]:
    ax = axes[0 if stat_name == 'R_fold' else 1]
    for idx, (name, label, color) in enumerate(configs):
        ckpt = outputs_dir / f'ablation_{name}.npz'
        if not ckpt.exists():
            continue
        data = np.load(ckpt)
        sim = data[stat_name]
        if len(sim) < 5:
            continue
        res = compare(obs_vals, sim, n_bootstrap=1000)
        ax.hist(sim, bins=30, range=(0, 1), density=True,
                histtype='step', linewidth=1.5, color=color, alpha=0.8)
        ax.hist(sim, bins=30, range=(0, 1), density=True,
                histtype='stepfilled', color=color, alpha=0.12)
        label_text = f'{label}  (KS p={res["ks_p_value"]:.3f})'
        # dummy for legend
        ax.plot([], [], color=color, linewidth=2, label=label_text)
    for r in obs_vals:
        ax.axvline(r, color='#d7191c', linewidth=1, alpha=0.5)
    ax.axvline(obs_vals[0], color='#d7191c', linewidth=1.5, alpha=0.8, label=f'Observed (N={len(obs_vals)})')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=7, loc='upper right')

fig.suptitle('Ablation: R_fold and R_cusp distributions across model complexity', fontsize=14, y=1.02)
fig.tight_layout()
out_path = outputs_dir / 'ablation_comparison.png'
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f'Saved {out_path}')
plt.close(fig)
