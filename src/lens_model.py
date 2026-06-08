import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver

COSMO_H0 = 67.4
COSMO_OMEGA_M = 0.315
COSMO_OMEGA_LAMBDA = 0.685
SPEED_OF_LIGHT = 299792.458


class MacroLens:
    def __init__(self):
        self._model = LensModel(['SIE', 'SHEAR'])
        self._solver = LensEquationSolver(self._model)
        self._kwargs_macro = []
        self._subhalo_list = []

    @staticmethod
    def _ang_dist_dc(z, z_ref=0):
        from scipy.integrate import quad
        E = lambda zp: np.sqrt(COSMO_OMEGA_M * (1 + zp) ** 3 + COSMO_OMEGA_LAMBDA)
        DH = 2997.92458 / COSMO_H0 * 100.0
        dc = DH * quad(lambda zp: 1.0 / E(zp), z_ref, z)[0]
        return dc

    def _einstein_radius(self, sigma_v, z_l, z_s):
        D_s = self._ang_dist_dc(z_s)
        D_d = self._ang_dist_dc(z_l)
        D_ds = self._ang_dist_dc(z_s, z_l)
        theta_E_rad = 4.0 * np.pi * (sigma_v / SPEED_OF_LIGHT) ** 2 * D_ds / D_s
        theta_E_arcsec = theta_E_rad * 206264.806247
        return theta_E_arcsec

    def build(self, z_l, z_s, sigma_v, q, phi_lens, gamma_ext, theta_gamma):
        theta_E = self._einstein_radius(sigma_v, z_l, z_s)
        e = (1.0 - q) / (1.0 + q)
        e1 = e * np.cos(2.0 * phi_lens)
        e2 = e * np.sin(2.0 * phi_lens)
        self._kwargs_macro = [
            {'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0.0, 'center_y': 0.0},
            {'gamma1': gamma_ext * np.cos(2.0 * theta_gamma),
             'gamma2': gamma_ext * np.sin(2.0 * theta_gamma)},
        ]
        self._rebuild_model()
        return theta_E

    def _rebuild_model(self):
        profiles = ['SIE', 'SHEAR']
        kwargs_all = list(self._kwargs_macro)
        for sub in self._subhalo_list:
            profiles.append('TNFW')
            kwargs_all.append(sub)
        self._model = LensModel(profiles)
        self._solver = LensEquationSolver(self._model)

    def add_substructure(self, tnfw_kwargs_list):
        self._subhalo_list.extend(tnfw_kwargs_list)
        self._rebuild_model()

    def clear_substructure(self):
        self._subhalo_list = []
        self._rebuild_model()

    def solve(self, beta_x, beta_y, num_random=20, search_window=5):
        theta_x, theta_y = self._solver.find_bright_image(
            beta_x, beta_y,
            self._kwargs_macro + self._subhalo_list,
            min_distance=0.01,
            search_window=search_window,
            num_random=num_random,
            initial_guess_cut=True,
        )
        if len(theta_x) < 4:
            return None
        kwargs_all = self._kwargs_macro + self._subhalo_list
        mu = self._model.magnification(theta_x, theta_y, kwargs_all)
        parity = np.array([1 if m > 0 else -1 for m in mu])
        sort_idx = np.lexsort((theta_y, theta_x))
        return {
            'theta_x': np.array(theta_x)[sort_idx],
            'theta_y': np.array(theta_y)[sort_idx],
            'mu': np.array(mu)[sort_idx],
            'parity': parity[sort_idx],
            'n_images': len(theta_x),
        }

    def solve_macro_only(self, beta_x, beta_y):
        fast_solver = LensEquationSolver(LensModel(['SIE', 'SHEAR']))
        theta_x, theta_y = fast_solver.find_bright_image(
            beta_x, beta_y, self._kwargs_macro,
            min_distance=0.01, search_window=5, num_random=5,
            initial_guess_cut=True,
        )
        if len(theta_x) < 4:
            return None
        mu = LensModel(['SIE', 'SHEAR']).magnification(theta_x, theta_y, self._kwargs_macro)
        parity = np.array([1 if m > 0 else -1 for m in mu])
        return {
            'theta_x': np.array(theta_x),
            'theta_y': np.array(theta_y),
            'mu': np.array(mu),
            'parity': parity,
            'n_images': len(theta_x),
        }
