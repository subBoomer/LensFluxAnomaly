import numpy as np
from data.radio_quads import build_observed_catalog
from src.statistic import compute_rfold

for sys in build_observed_catalog():
    parity = np.array(sys['parity'])
    x = np.array(sys['x_arcsec'])
    y = np.array(sys['y_arcsec'])
    f = np.array(sys['fluxes_mjy'])
    r = compute_rfold(x, y, parity, f)
    pos_idx = [i for i, p in enumerate(parity) if p > 0]
    r_dist = np.sqrt(x[pos_idx]**2 + y[pos_idx]**2)
    names = sys['image_ids']
    pos_names = [names[i] for i in pos_idx]
    print(f'{sys["name"]:>12s}: R_fold={r:.4f}  parity={list(parity)}  pos_pair={pos_names} radii={np.round(r_dist,3)}')
