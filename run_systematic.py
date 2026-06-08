import sys, numpy as np, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from data.radio_quads import build_observed_catalog
from src.statistic import compute_rfold
from src.comparison import compare

outputs_dir = Path(__file__).parent / 'outputs'

systems = build_observed_catalog()
R_obs = []
for s in systems:
    x = np.array(s['x_arcsec']); y = np.array(s['y_arcsec'])
    p = np.array(s['parity']); f = np.array(s['fluxes_mjy'])
    r = compute_rfold(x, y, p, f)
    if r >= 0: R_obs.append(r)
R_obs = np.array(R_obs)

print('=== SNR sensitivity ===')
print(f'{"SNR":>5} {"N_sim":>6} {"KS_p":>9} {"Wass":>9}')
print('-' * 33)

for snr in [20, 50, 100]:
    ckpt = outputs_dir / f'systest_snr_{snr}.npz'
    if ckpt.exists():
        R_sim = np.load(ckpt)['R_sim']
    else:
        r = subprocess.run([
            sys.executable, 'run_inference.py',
            '--n-systems', '500', '--force', '--parallel', '--snr', str(snr)
        ], capture_output=True, text=True, timeout=300)
        print(r.stdout.strip().split('\n')[-4:])
        R_sim = np.load(outputs_dir / 'inference_result.npz')['R_sim']
        np.savez(ckpt, R_sim=R_sim)
    res = compare(R_obs, R_sim, n_bootstrap=1000)
    print(f'{snr:>5d} {len(R_sim):>6d} {res["ks_p_value"]:>9.4f} {res["wasserstein_distance"]:>9.4f}')

print()

# Optional: concentration model test
print('=== Concentration model ===')
for model in ['duffy08', 'ishiyama21', 'diemer15']:
    ckpt = outputs_dir / f'systest_conc_{model}.npz'
    if ckpt.exists():
        R_sim = np.load(ckpt)['R_sim']
    else:
        # modify substructure.py CONCENTRATION_MODEL
        mod_path = Path(__file__).parent / 'src' / 'substructure.py'
        orig = mod_path.read_text(encoding='utf-8')
        mod = orig.replace("CONCENTRATION_MODEL = 'duffy08'", f"CONCENTRATION_MODEL = '{model}'")
        mod_path.write_text(mod, encoding='utf-8')
        try:
            r = subprocess.run([
                sys.executable, 'run_inference.py',
                '--n-systems', '500', '--force', '--parallel'
            ], capture_output=True, text=True, timeout=300)
            print(r.stdout.strip().split('\n')[-4:])
            R_sim = np.load(outputs_dir / 'inference_result.npz')['R_sim']
            np.savez(ckpt, R_sim=R_sim)
        finally:
            mod_path.write_text(orig, encoding='utf-8')
    res = compare(R_obs, R_sim, n_bootstrap=1000)
    print(f'{model:20s} {len(R_sim):>6d} {res["ks_p_value"]:>9.4f} {res["wasserstein_distance"]:>9.4f}')
