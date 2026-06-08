import sys, warnings, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.resolve()))
warnings.filterwarnings('ignore', category=UserWarning)

from src.statistic import compute_rfold
from src.perturbation_model import SimplePerturbationModel

N_SIM = 50000
SEED = 42


def main():
    print('=' * 60)
    print('Selection Bias Test')
    print('=' * 60)
    print('CLASS selection criteria:')
    print('  * Flux limit: 30 mJy (parent survey)')
    print('  * Resolution: 0.2" VLA detected multi-component')
    print('  * Spectral index: flat (alpha > -0.5) - source property')
    print('  * No flux-ratio selection')
    print()
    print('Key question: does CLASS selection correlate with R_fold?')
    print('If selection picks specific R_fold values, the catalog is biased.')
    print()

    model = SimplePerturbationModel()
    rng = np.random.default_rng(SEED)

    rfold_list = []
    selected_list = []
    flux_tot_list = []
    min_sep_list = []

    for i in range(N_SIM):
        result = model.realise(seed=rng.integers(0, 2**31), return_all=True)
        if result is None or result[0] is None:
            continue

        rmin_val, x, y, F, *_ = result
        parity = np.array([1, -1, 1, -1])
        rf = compute_rfold(x, y, parity.astype(float), F)
        if rf < 0:
            continue
        rfold_list.append(rf)

        # Approximate lensed flux: F are normalized (sum=1).
        # Typical CLASS quad: source flux ~10 mJy, magnification ~5 => total ~50 mJy.
        source_flux = 10.0
        macro_mu_tot = 5.0
        total_flux = macro_mu_tot * source_flux  # ~50 mJy

        min_sep = float('inf')
        n = len(x)
        for a in range(n):
            for b in range(a + 1, n):
                d = np.hypot(x[a] - x[b], y[a] - y[b])
                min_sep = min(min_sep, d)
        min_sep_list.append(min_sep)
        flux_tot_list.append(total_flux)

        # CLASS selection criteria (simplified):
        # Flux cut: total > 30 mJy
        # Resolution: min_sep > 0.15" (slightly below VLA 0.2" beam for detectability)
        passes = (total_flux > 30.0) and (min_sep > 0.15)
        selected_list.append(passes)

        if (i + 1) % 10000 == 0:
            print(f'  {i+1}/{N_SIM}', flush=True)

    rfold_arr = np.array(rfold_list)
    selected_arr = np.array(selected_list)
    min_sep_arr = np.array(min_sep_list)

    print(f'\nTotal CDM realizations: {len(rfold_arr)}')
    print(f'Pass CLASS selection:   {selected_arr.sum()} ({100*selected_arr.mean():.1f}%)')
    print(f'Fail (flux < 30 mJy):   {0} (all pass flux cut at nominal flux)')
    print(f'Fail (sep < 0.15\"):      {(min_sep_arr < 0.15).sum()} ({100*(min_sep_arr < 0.15).mean():.1f}%)')

    sel = rfold_arr[selected_arr]
    unsel = rfold_arr[~selected_arr]

    print(f'\n{"Metric":25s} {"All":>10s} {"Selected":>10s} {"Unselected":>10s}')
    print('-' * 57)
    print(f'{"N":25s} {len(rfold_arr):>10d} {len(sel):>10d} {len(unsel):>10d}')
    print(f'{"Mean R_fold":25s} {np.mean(rfold_arr):>10.4f} {np.mean(sel):>10.4f} {np.mean(unsel):>10.4f}')
    print(f'{"Std R_fold":25s} {np.std(rfold_arr):>10.4f} {np.std(sel):>10.4f} {np.std(unsel):>10.4f}')
    print(f'{"P(R>0.2)":25s} {np.mean(rfold_arr>0.2):>10.4f} {np.mean(sel>0.2):>10.4f} {np.mean(unsel>0.2):>10.4f}')

    delta = np.mean(sel) - np.mean(rfold_arr)
    obs_cdm_gap = 0.436 - 0.165
    print(f'\nDelta mean (selected - all): {delta:+.4f}')
    print(f'Obs-CDM gap (R_fold): {obs_cdm_gap:.3f}')
    print(f'Selection bias / anomaly: {abs(delta)/obs_cdm_gap*100:.1f}%')

    print(f'\n--- Interpretation ---')
    if abs(delta) < 0.02:
        print('Selection bias < 0.02 in mean R_fold → NEGLIGIBLE.')
        print('CLASS selection does not correlate with R_fold.')
    else:
        print(f'Selection bias = {delta:.3f} in mean R_fold.')
        print(f'Direction: {"selected R_fold > all" if delta > 0 else "selected R_fold < all"}')
        print(f'Bias fraction: {abs(delta)/obs_cdm_gap*100:.1f}% of anomaly signal.')
        print('Needs to be accounted for in inference.')

    np.savez('outputs/selection_bias_test.npz',
             rfold_all=rfold_arr, rfold_selected=sel, min_sep=min_sep_arr)

    print(f'\nSaved: outputs/selection_bias_test.npz')


if __name__ == '__main__':
    main()
