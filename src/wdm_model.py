import numpy as np
from src.perturbation_model import SimplePerturbationModel


class WDMPerturbationModel(SimplePerturbationModel):
    """Simple perturbation model with WDM mass function suppression.

    Suppresses the subhalo mass function below the half-mode mass:
        n_WDM(m) / n_CDM(m) = 1 / (1 + (m_hm / m)^beta)

    Half-mode mass for thermal WDM particle mass m_x (keV):
        m_hm = 1e8 * (m_x / 3 keV)^{-3.33} M_sun
    """

    def __init__(self, m_wdm_keV=3.0, beta=2.5, **kwargs):
        super().__init__(**kwargs)
        self.m_wdm_keV = m_wdm_keV
        self.beta = beta
        self.m_hm = 1e8 * (m_wdm_keV / 3.0) ** (-3.33)

    def suppression_factor(self, masses):
        return 1.0 / (1.0 + (self.m_hm / masses) ** self.beta)

    def draw_subhalos(self, rng):
        n = rng.poisson(self.n_sub_mean)
        if n == 0:
            return np.array([]), np.array([]), np.array([])
        m_min_val = self.m_min ** (1 - self.alpha)
        m_max_val = self.m_max ** (1 - self.alpha)
        masses = (rng.uniform(size=n) * (m_max_val - m_min_val) + m_min_val) ** (1 / (1 - self.alpha))
        keep = rng.uniform(size=n) < self.suppression_factor(masses)
        masses = masses[keep]
        if len(masses) == 0:
            return np.array([]), np.array([]), np.array([])
        r_sub = rng.exponential(scale=1.5, size=len(masses))
        theta_sub = rng.uniform(0, 2 * np.pi, size=len(masses))
        x_sub = r_sub * np.cos(theta_sub)
        y_sub = r_sub * np.sin(theta_sub)
        return masses, x_sub, y_sub
