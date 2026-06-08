import warnings
import numpy as np
from colossus.cosmology import cosmology
from colossus.halo.concentration import concentration
from colossus.lss import mass_function

cos = cosmology.setCosmology('planck18')
warnings.filterwarnings('ignore', category=UserWarning,
                        module='colossus.halo.concentration')

F_SUB_DMO = 0.005
MASS_SLOPE = -1.9
M_MIN = 1e6
M_MAX = 1e10
LAMBDA_DMO = 50.0
CONCENTRATION_MODEL = 'duffy08'


def _rho_crit(z):
    Ez = np.sqrt(cos.Om0 * (1 + z) ** 3 + cos.Ode0)
    rho_crit0 = 2.77536627e11 * cos.h ** 2
    return rho_crit0 * Ez ** 2


class SubhaloPopulation:
    def __init__(self, f_sub=None):
        self.f_sub = F_SUB_DMO if f_sub is None else f_sub
        self.lambda_dmo = LAMBDA_DMO

    def _sample_mass(self, rng):
        slope = MASS_SLOPE + 1.0
        m_low = M_MIN ** slope
        m_high = M_MAX ** slope
        u = rng.uniform(0, 1)
        m = (m_low + u * (m_high - m_low)) ** (1.0 / slope)
        return m

    def realise(self, theta_E, z_l, rng):
        lambda_eff = self.lambda_dmo * self.f_sub / F_SUB_DMO
        n_sub = rng.poisson(lambda_eff)
        subhalos = []
        r_max = theta_E
        h = cos.h
        for _ in range(n_sub):
            m = self._sample_mass(rng)
            m_h = m / h
            r = rng.uniform(0, r_max)
            angle = rng.uniform(0, 2 * np.pi)
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            try:
                c = concentration(m_h, '200c', z_l, model=CONCENTRATION_MODEL)
            except Exception:
                c = 10.0
            rho_c = _rho_crit(z_l)
            r_delta = (3.0 * m / (4.0 * np.pi * 200.0 * rho_c)) ** (1.0 / 3.0)
            rs = r_delta / c
            r_trunc = max(rs * 0.5, r_delta / 10.0)
            subhalos.append({
                'alpha_Rs': rs,
                'Rs': rs,
                'r_trunc': r_trunc,
                'center_x': x,
                'center_y': y,
            })
        return subhalos


class LOSPopulation:
    def __init__(self):
        pass

    def realise(self, z_l, z_s, rng):
        D_C_s = cos.comovingDistance(0, z_s)
        D_C_ref = cos.comovingDistance(0, 2.0)
        z_factor = (1.0 + z_l) / 1.5
        sigma_kappa = 0.030 * np.sqrt(D_C_s / D_C_ref) * z_factor

        sigma_gamma = sigma_kappa * 0.8
        corr = 0.4
        cov = np.array([
            [sigma_kappa ** 2, corr * sigma_kappa * sigma_gamma, 0.0],
            [corr * sigma_kappa * sigma_gamma, sigma_gamma ** 2, 0.0],
            [0.0, 0.0, sigma_gamma ** 2],
        ])
        kappa, gamma_1, gamma_2 = rng.multivariate_normal(np.zeros(3), cov)
        return {'kappa': kappa, 'gamma_1': gamma_1, 'gamma_2': gamma_2}
