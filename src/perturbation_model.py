import numpy as np


class SimplePerturbationModel:
    """Minimal subhalo perturbation model — no ray tracing.

    Generates abstract quad image positions for an SIE-like lens, then
    perturbs the macro fluxes using a phenomenological kernel:

        δ_ik = ε · m_k / (r_ik² + r_c²)

    The macro (unperturbed) flux is set to equal F_i⁰ = 0.25 for all 4 images.
    Subhalo masses are drawn from dN/dm ∝ m^{-alpha}, placed in 2D projection
    following an NFW-tracing spatial distribution.
    """

    def __init__(
        self,
        epsilon=5e-9,
        r_c_kpc=1.0,
        m_min=1e6,
        m_max=1e9,
        alpha=1.9,
        sigma_noise=0.05,
        n_sub_mean=200,
    ):
        self.epsilon = epsilon
        self.r_c_kpc = r_c_kpc
        self.m_min = m_min
        self.m_max = m_max
        self.alpha = alpha
        self.sigma_noise = sigma_noise
        self.n_sub_mean = n_sub_mean
        self._rng = None

    def _ensure_rng(self, seed):
        if self._rng is None or seed is not None:
            self._rng = np.random.default_rng(seed)
        return self._rng

    def generate_abstract_quad(self, rng):
        """Generate 4 symmetric image positions for an abstract SIE-like lens.

        Returns
        -------
        x, y : ndarrays of shape (4,)
            Image positions in arcsec.
        """
        theta_E = rng.lognormal(mean=0.0, sigma=0.4) * 1.0
        theta_E = max(theta_E, 0.3)
        theta_E = min(theta_E, 2.5)
        e = rng.uniform(0.1, 0.5)
        q = (1.0 - e) / (1.0 + e)
        phi = rng.uniform(0, np.pi)
        r1 = theta_E * (1.0 - e * 0.5)
        r2 = theta_E * (1.0 + e * 0.3)
        angles = rng.uniform(0, np.pi / 2) + np.array([0, np.pi / 2, np.pi, 3 * np.pi / 2])
        radii = np.array([r1, r2, r1, r2])
        scatter = rng.uniform(0.85, 1.15, size=4)
        x = radii * scatter * np.cos(angles + phi)
        y = radii * scatter * np.sin(angles + phi)
        return x, y

    def draw_subhalos(self, rng):
        """Draw subhalo masses and 2D positions.

        Returns
        -------
        masses : ndarray of shape (N,)
            Subhalo masses in M_sun.
        x_sub, y_sub : ndarrays of shape (N,)
            Subhalo positions in arcsec (centered on lens).
        """
        n = rng.poisson(self.n_sub_mean)
        if n == 0:
            return np.array([]), np.array([]), np.array([])
        m_min_val = self.m_min ** (1 - self.alpha)
        m_max_val = self.m_max ** (1 - self.alpha)
        masses = (rng.uniform(size=n) * (m_max_val - m_min_val) + m_min_val) ** (1 / (1 - self.alpha))
        r_sub = rng.exponential(scale=1.5, size=n)
        theta_sub = rng.uniform(0, 2 * np.pi, size=n)
        x_sub = r_sub * np.cos(theta_sub)
        y_sub = r_sub * np.sin(theta_sub)
        return masses, x_sub, y_sub

    def _compute_perturbations(self, x_img, y_img, masses, x_sub, y_sub):
        """Compute δ_i = Σ_k ε · m_k / (r_ik² + r_c²) for each image."""
        r_c_arcsec = self.r_c_kpc * 0.25
        deltas = np.zeros(4)
        if len(masses) == 0:
            return deltas
        for i in range(4):
            r2 = (x_img[i] - x_sub) ** 2 + (y_img[i] - y_sub) ** 2
            deltas[i] = np.sum(self.epsilon * masses / (r2 + r_c_arcsec ** 2))
        return deltas

    def realise(self, seed=None, return_all=False):
        """Generate one full realisation: quad positions → perturbed fluxes → R_min.

        Parameters
        ----------
        seed : int or None
        return_all : bool
            If True, return (rmin, x_img, y_img, F_perturbed, masses, x_sub, y_sub).

        Returns
        -------
        rmin : float or None
        """
        rng = self._ensure_rng(seed)
        x_img, y_img = self.generate_abstract_quad(rng)
        masses, x_sub, y_sub = self.draw_subhalos(rng)
        deltas = self._compute_perturbations(x_img, y_img, masses, x_sub, y_sub)
        F_unpert = np.ones(4) * 0.25
        F_pert = F_unpert * (1.0 + deltas)
        F_pert = np.maximum(F_pert, 0.0)
        noise = rng.normal(0, self.sigma_noise, size=4)
        F_obs = np.maximum(F_pert * (1.0 + noise), 0.0)
        F_obs /= np.sum(F_obs)
        x_c, y_c = np.mean(x_img), np.mean(y_img)
        r = np.hypot(x_img - x_c, y_img - y_c)
        deltas_r = abs(r[:, None] - r[None, :])
        best = None
        for i in range(4):
            for j in range(i + 1, 4):
                if deltas_r[i, j] < 0.2:
                    val = abs(F_obs[i] - F_obs[j]) / (F_obs[i] + F_obs[j])
                    if best is None or val < best:
                        best = val
        if return_all:
            return best, x_img, y_img, F_obs, masses, x_sub, y_sub
        return best

    def sample(self, n_realisations, seed=42, verbose=True):
        """Generate many realisations and return array of R_min values.

        Parameters
        ----------
        n_realisations : int
        seed : int
        verbose : bool

        Returns
        -------
        rmin_values : ndarray
            Array of R_min values (None entries filtered).
        """
        rng = self._ensure_rng(seed)
        results = []
        for i in range(n_realisations):
            rmin = self.realise(seed=rng.integers(0, 2 ** 31))
            if rmin is not None:
                results.append(rmin)
            if verbose and (i + 1) % 5000 == 0:
                print(f'  Sim {i + 1}/{n_realisations} ({len(results)} valid)')
        return np.array(results)
