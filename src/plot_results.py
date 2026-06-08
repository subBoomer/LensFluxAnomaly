import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

outputs_dir = Path(__file__).parent.parent / 'outputs'
data = np.load(outputs_dir / 'inference_result.npz')

R_obs = data['R_obs']
R_sim = data['R_sim']
ks_p = float(data['ks_p_value'][0])
wass = float(data['wasserstein_distance'][0])

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(R_sim, bins=40, range=(0, 1), density=True,
        alpha=0.7, color='#2c7bb6', label=f'Simulated (N={len(R_sim)})')
for r in R_obs:
    ax.axvline(r, color='#d7191c', linewidth=1.5, alpha=0.8)
ax.axvline(R_obs[0], color='#d7191c', linewidth=1.5, alpha=0.8, label=f'Observed (N={len(R_obs)})')

ax.set_xlabel(r'$R_{\mathrm{fold}}$', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title(rf'$R_{{\mathrm{{fold}}}}$ Distribution: KS $p={ks_p:.3f}$, W={wass:.3f}', fontsize=13)
ax.legend(fontsize=11)
ax.set_xlim(0, 1)
fig.tight_layout()

plot_path = outputs_dir / 'rfold_distribution.png'
fig.savefig(plot_path, dpi=150)
print(f'Saved plot to {plot_path}')
plt.close(fig)
