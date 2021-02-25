import numpy as np
import scipy.optimize as opt
from math import pi, sqrt, factorial
from scipy.integrate import simps
from scipy.special import sph_harm
from numba import njit, prange, int32, float64, complex128


def sphere_integration(dphi, theta, func):
    """
    Integrate `func` in spherical coordinates
    """
    func_theta = np.array(
        [phi_integ for phi_integ in simps(func, dx=dphi, axis=1)]
    )
    return simps(np.sin(theta) * func_theta, dx=theta[1] - theta[0])


@njit(
    (int32, int32, int32, complex128[:], complex128[:, :, :], complex128[:, :])
)
def set_state(ntht, nphi, n_sph_harm, coef, sph_harm_arr, state):
    """
    Auxiliar function to compute state from coefficients
    in spherical harmonics basis expansion

    output
    ------
    `state` : ``2D numpy.array``

    """
    cusum = complex(0, 0)
    for i in prange(ntht):
        for j in prange(nphi):
            cusum = 0.0 + 0.0j
            for k in prange(n_sph_harm):
                cusum = cusum + coef[k] * sph_harm_arr[k, i, j]
            state[i, j] = cusum


class StateSphericalBasis:
    """
    Class to provide interface in using constraint minimization in spherical
    harmonic basis with scipy.optimize.minimize and 'SLSQP' method. In class
    initialization, minimum and maximum values for 'l' angular momentum must
    be provided, as well as the grid points in the sphere

    Parameters
    ----------
    `lmax` : ``int``
        maximum value for l angular momentum in basis expansion
    `lmin` : ``int``
        minimum value for l angular momentum in basis expansion <= `lmax`
    `m_min` : ``int``
        minimum ABSOLUTE value for azimuthal number <= `lmin`
    `n_phi_pts` : ``int`` (default 101)
        number of grid points to discretize theta angle
    `n_tht_pts` : ``int`` (default 101)
        number of grid points to discretize phi angle

    Usual parameters in methods
    `frac_a`, `frac_b`, `ga`, `gb`, `gab`

    """

    def __init__(
        self, lmax, lmin=0, m_min=0, n_phi_pts=101, n_tht_pts=101, seed=1000
    ):
        self.lmax = lmax
        self.lmin = lmin
        self.m_min = m_min
        self.phi_pts = np.linspace(0, 2 * pi, n_phi_pts)
        self.tht_pts = np.linspace(0, pi, n_tht_pts)
        self.phi_grid, self.tht_grid = np.meshgrid(self.phi_pts, self.tht_pts)
        self.n_sph_harm = (lmax + 1 - m_min + lmin + 1 - m_min) * (
            lmax - lmin + 1
        )
        if m_min == 0:
            self.n_sph_harm = self.n_sph_harm - (lmax - lmin + 1)
        self.sph_harm_arr = np.empty(
            [self.n_sph_harm, n_tht_pts, n_phi_pts], dtype=np.complex128
        )
        i = 0
        for l in range(lmin, lmax + 1):
            default_m = np.arange(-l, l + 1)
            consider_m = default_m[abs(default_m) >= m_min]
            for m in consider_m:
                self.sph_harm_arr[i] = sph_harm(
                    m, l, self.phi_grid, self.tht_grid
                )
                i = i + 1
        self.seed = seed

    def __kinect_energy(self, coef_a, coef_b, frac_a, frac_b):
        k = 0
        kin = 0.0
        abs2_coef_a = abs(coef_a) ** 2
        abs2_coef_b = abs(coef_b) ** 2
        for l in range(self.lmin, self.lmax + 1):
            default_m = np.arange(-l, l + 1)
            consider_m = default_m[abs(default_m) >= self.m_min]
            for m in consider_m:
                kin = kin + 0.5 * l * (l + 1) * (
                    frac_a * abs2_coef_a[k] + frac_b * abs2_coef_b[k]
                )
                k = k + 1
        return kin

    def __energy_functional(self, c, frac_a, frac_b, ga, gb, gab):
        coef_a = (
            c[: self.n_sph_harm]
            + 1.0j * c[self.n_sph_harm : 2 * self.n_sph_harm]
        )
        coef_b = (
            c[2 * self.n_sph_harm : 3 * self.n_sph_harm]
            + 1.0j * c[3 * self.n_sph_harm :]
        )
        state_a = np.empty(self.sph_harm_arr[0].shape, dtype=np.complex128)
        state_b = np.empty(self.sph_harm_arr[0].shape, dtype=np.complex128)
        set_state(
            self.tht_pts.size,
            self.phi_pts.size,
            self.n_sph_harm,
            coef_a,
            self.sph_harm_arr,
            state_a,
        )
        set_state(
            self.tht_pts.size,
            self.phi_pts.size,
            self.n_sph_harm,
            coef_b,
            self.sph_harm_arr,
            state_b,
        )
        inter_pot = (
            0.5 * ga * frac_a * abs(state_a) ** 4
            + 0.5 * gb * frac_b * abs(state_b) ** 4
            + gab
            * sqrt(frac_a * frac_b)
            * abs(state_a) ** 2
            * abs(state_b) ** 2
        )
        inter_energy = sphere_integration(
            self.phi_pts[1] - self.phi_pts[0], self.tht_pts, inter_pot
        )
        kin = self.__kinect_energy(coef_a, coef_b, frac_a, frac_b)
        return kin + inter_energy

    def __norm_constraint_a(self, c):
        return (c[: int(c.size / 2)] ** 2).sum() - 1

    def __norm_constraint_b(self, c):
        return (c[int(c.size / 2) :] ** 2).sum() - 1

    def __separate_coef(self, c):
        ca = (
            c[: self.n_sph_harm]
            + 1.0j * c[self.n_sph_harm : 2 * self.n_sph_harm]
        )
        cb = (
            c[2 * self.n_sph_harm : 3 * self.n_sph_harm]
            + 1.0j * c[3 * self.n_sph_harm :]
        )
        return ca, cb

    def __state(self, coef):
        psi = np.zeros(
            [self.tht_pts.size, self.phi_pts.size], dtype=np.complex128
        )
        i = 0
        for l in range(self.lmin, self.lmax + 1):
            default_m = np.arange(-l, l + 1)
            consider_m = default_m[abs(default_m) >= self.m_min]
            for m in consider_m:
                psi = psi + coef[i] * self.sph_harm_arr[i]
                i = i + 1
        return psi

    def minimize_energy(self, frac_a, frac_b, ga, gb, gab, err_tol=1e-4):
        """
        Minimize Gross-Pitaevskii energy functional in the domain defined
        in class initialization using spherical harmonics basis expansion

        Parameters
        ----------
        `frac_a` : ``float64``
            fraction of number of atoms of species A
        `frac_b` : ``float64``
            fraction of number of atoms of species B
        `ga` : ``float64``
            interaction parameter of species A
        `gb` : ``float64``
            interaction parameter of species B
        `gab` : ``float64``
            interspecies interaction parameter
        `err_tol` : ``float64`` default 1E-4
            error allowed in optimization method (optional)

        Return
        ------
        ``tuple`` : (``float``, ``2D numpy.array``, ``2D numpy.array``)
            (energy, state_a, state_b)

        """
        c0 = np.zeros(2 * self.n_sph_harm, dtype=np.complex128)
        rand_generator = np.random.default_rng(seed=self.seed)
        w_real = rand_generator.random(2 * self.n_sph_harm) - 0.5
        w_imag = rand_generator.random(2 * self.n_sph_harm) - 0.5
        i = 0
        for l in range(self.lmin, self.lmax + 1):
            default_m = np.arange(-l, l + 1)
            consider_m = default_m[abs(default_m) >= self.m_min]
            for m in consider_m:
                c0[i] = (w_real[i] + 1.0j * w_imag[i]) / sqrt(factorial(l + 1))
                j = i + self.n_sph_harm
                c0[j] = (w_real[j] + 1.0j * w_imag[j]) / sqrt(factorial(l + 1))
                i = i + 1
        c0_spec_a = c0[: self.n_sph_harm]
        c0_spec_b = c0[self.n_sph_harm :]
        c0_spec_a = c0_spec_a / sqrt((abs(c0_spec_a) ** 2).sum())
        c0_spec_b = c0_spec_b / sqrt((abs(c0_spec_b) ** 2).sum())
        c0_reshape = np.concatenate(
            [c0_spec_a.real, c0_spec_a.imag, c0_spec_b.real, c0_spec_b.imag]
        )
        cons = [
            {"type": "eq", "fun": self.__norm_constraint_a},
            {"type": "eq", "fun": self.__norm_constraint_b},
        ]
        extra = {"ftol": err_tol, "disp": False, "maxiter": 5000}
        res = opt.minimize(
            self.__energy_functional,
            c0_reshape,
            method="SLSQP",
            args=(frac_a, frac_b, ga, gb, gab),
            constraints=cons,
            options=extra,
        )
        ca, cb = self.__separate_coef(res.x)
        return res.fun, self.__state(ca), self.__state(cb)

    def sweep_interspecies(self, frac_a, frac_b, ga, gb, gab_list):
        """
        Compute states using minimization routine for a list of
        interspecies interaction strength `gab_list` with other
        parameters fixed
        """
        results = []
        for gab in gab_list:
            res = self.minimize_energy(frac_a, frac_b, ga, gb, gab)
            results.append(res)
        return results
