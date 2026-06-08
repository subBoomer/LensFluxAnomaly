import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

COSMO_H0 = 67.4
COSMO_OMEGA_M = 0.315
COSMO_OMEGA_LAMBDA = 0.685


class LensPopulation:
    def __init__(self, rng: np.random.Generator | None = None):
        self._rng = rng
        self._z_grid = None
        self._dVdz_grid = None
        self._cdf_interp = None
        self._build_redshift_cdf()

    def _build_redshift_cdf(self):
        z_grid = np.linspace(0.01, 1.5, 500)
        dVdz = np.array([self._dVdz(z) for z in z_grid])
        pdf = dVdz * (1.0 + z_grid) ** 2
        pdf /= np.trapezoid(pdf, z_grid)
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        self._z_grid = z_grid
        self._dVdz_grid = dVdz
        self._cdf_interp = interp1d(cdf, z_grid, bounds_error=False, fill_value=(z_grid[0], z_grid[-1]))

    @staticmethod
    def _E(z):
        return np.sqrt(COSMO_OMEGA_M * (1 + z) ** 3 + COSMO_OMEGA_LAMBDA)

    @staticmethod
    def _dVdz(z):
        DH = 2997.92458 / COSMO_H0 * 100.0
        Ez = LensPopulation._E(z)
        dc = DH * quad(lambda zp: 1.0 / LensPopulation._E(zp), 0, z)[0]
        return DH * (1 + z) ** 2 * dc ** 2 / Ez

    def _sample_redshift(self):
        u = self._rng.uniform(0, 1)
        return float(self._cdf_interp(u))

    @staticmethod
    def _axis_ratio(rng):
        q = rng.normal(0.8, 0.15)
        return np.clip(q, 0.3, 1.0)

    def sample(self, rng: np.random.Generator | None = None):
        if rng is not None:
            self._rng = rng
        rng = self._rng
        z_l = self._sample_redshift()
        z_s = rng.uniform(1.0, 3.0)
        log_sigma_v = rng.normal(np.log(250.0), 0.15)
        sigma_v = np.exp(log_sigma_v)
        q = self._axis_ratio(rng)
        gamma_ext = rng.exponential(0.05)
        theta_gamma = rng.uniform(0, np.pi)
        source_flux = rng.lognormal(-1.0, 0.5)
        phi_lens = rng.uniform(0, np.pi)
        return {
            'z_l': z_l,
            'z_s': z_s,
            'sigma_v': sigma_v,
            'q': q,
            'phi_lens': phi_lens,
            'gamma_ext': gamma_ext,
            'theta_gamma': theta_gamma,
            'source_flux': source_flux,
            'beta_x': None,
            'beta_y': None,
        }
