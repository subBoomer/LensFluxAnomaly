import numpy as np
d = np.load('outputs/checkpoint.npz')
print(f'completed={d["completed"]}, n_R_sim={len(d["R_sim"])}')
