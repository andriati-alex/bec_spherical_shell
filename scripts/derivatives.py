import numpy as np
from scipy.integrate import simpson
from scipy.special import gamma
from scipy.optimize import minimize_scalar, OptimizeResult
from typing import Callable, Any


def der_reflect(f: np.ndarray, dx: float, azi: int) -> np.ndarray:
    n = f.size
    r = 1.0 / (12 * dx)
    sign = 1 - 2 * (abs(azi) % 2)
    dfdx = np.empty(n, dtype=f.dtype)
    dfdx[2 : n - 2] = r * (
        f[: n - 4] - f[4:] + 8 * (f[3 : n - 1] - f[1 : n - 3])
    )
    # set the boundary values
    dfdx[0] = r * ((sign - 1) * f[2] + 8 * (1 - sign) * f[1])
    dfdx[1] = r * (sign * f[1] - f[3] + 8 * (f[2] - f[0]))
    dfdx[n - 2] = r * (f[n - 4] - sign * f[n - 2] + 8 * (f[n - 1] - f[n - 3]))
    dfdx[n - 1] = r * ((1 - sign) * f[n - 3] + 8 * (sign - 1) * f[n - 2])
    return dfdx


def der2_reflect(f: np.ndarray, dx: float, azi: int) -> np.ndarray:
    n = f.size
    r = 1.0 / (12 * dx * dx)
    sign = 1 - 2 * (abs(azi) % 2)
    d2fdx = np.empty(n, dtype=f.dtype)
    d2fdx[2 : n - 2] = r * (
        -(f[: n - 4] + f[4:])
        + 16 * (f[3 : n - 1] + f[1 : n - 3])
        - 30 * f[2 : n - 2]
    )
    # set the boundary values
    d2fdx[0] = r * (-(1 + sign) * f[2] + 16 * (1 + sign) * f[1] - 30 * f[0])
    d2fdx[1] = r * (-(f[3] + sign * f[1]) + 16 * (f[2] + f[0]) - 30 * f[1])
    d2fdx[n - 1] = r * (
        -(1 + sign) * f[n - 3] + 16 * (1 + sign) * f[n - 2] - 30 * f[n - 1]
    )
    d2fdx[n - 2] = r * (
        -(f[n - 4] + sign * f[n - 2])
        + 16 * (f[n - 1] + f[n - 3])
        - 30 * f[n - 2]
    )
    return d2fdx


def energy_func_tht(
    tht: np.ndarray, f: np.ndarray, g: float, azi: int = 0
) -> float:
    n = tht.size
    dtht = tht[1] - tht[0]
    dfdx = der_reflect(f, dtht, azi)
    regular_azi_term = np.zeros(n, dtype=f.dtype)
    regular_azi_term[1 : n - 1] = (
        abs(azi * f[1 : n - 1] / np.sin(tht[1 : n - 1])) ** 2
    )
    # 0.5 comes from p^2 / (2 m) in dimensionless units
    kinect_den = 0.5 * abs(dfdx) ** 2 + 0.5 * regular_azi_term
    return simpson(
        np.sin(tht) * (kinect_den + 0.5 * g * abs(f) ** 4),
        dx=tht[1] - tht[0],
    )


def mu_func_tht(
    tht: np.ndarray, f: np.ndarray, g: float, azi: int = 0
) -> float:
    n = tht.size
    dtht = tht[1] - tht[0]
    dfdx = der_reflect(f, dtht, azi)
    regular_azi_term = np.zeros(n, dtype=f.dtype)
    regular_azi_term[1 : n - 1] = (
        abs(azi * f[1 : n - 1] / np.sin(tht[1 : n - 1])) ** 2
    )
    # 0.5 comes from p^2 / (2 m) in dimensionless units
    kinect_den = 0.5 * abs(dfdx) ** 2 + 0.5 * regular_azi_term
    return simpson(
        np.sin(tht) * (kinect_den + g * abs(f) ** 4),
        dx=tht[1] - tht[0],
    )


def gp_residue_tht(
    mu: float, tht: np.ndarray, f: np.ndarray, g: float, azi: int = 0
) -> float:
    n = f.size
    dtht = tht[1] - tht[0]
    dfdx = der_reflect(f, dtht, azi)
    d2fdx = der2_reflect(f, dtht, azi)
    regular_azi_term = np.zeros(n, dtype=f.dtype)
    regular_azi_term[1 : n - 1] = (
        -(azi ** 2) * f[1 : n - 1] / np.sin(tht[1 : n - 1]) ** 2
    )
    regular_der1 = np.zeros(n, dtype=f.dtype)
    regular_der1[1 : n - 1] = (
        dfdx[1 : n - 1] * np.cos(tht[1 : n - 1]) / np.sin(tht[1 : n - 1])
    )
    # 0.5 comes from p^2 / (2 m) in dimensionless units
    gp_f = (
        -0.5 * (d2fdx + regular_der1 + regular_azi_term)
        + g * (abs(f) ** 2) * f
    )
    return (
        np.sqrt(
            simpson(np.sin(tht) * abs(mu * f - gp_f) ** 2, dx=tht[1] - tht[0])
        )
        / mu
    )


def general_energy_min(
    func: Callable[..., np.ndarray],
    energy_kwargs: dict[str, Any],
    bounds: tuple[int, int] = (0, 100),
    func_args: tuple = (),
) -> OptimizeResult:
    """Wrapper to use callable for energy minimization (1 parameter)

    Use the function `scipy.optimize.minimize_scalar` to minimize the
    Goss-Pitaevskii energy functional. The main concern is to provide
    a function that is normalized regardless of the variational param
    chosen

    Parameters:
    -----------
    `func`: client function to generate a normalized state on the grid
    `e_kwargs`: all arguments of `energy_func_tht` but `f`
    `bounds`: Two numbers to restrict the minimization interval
    `func_args`: additional arguments required in `func`

    Return
    ------
    ``scipy.optimize.OptimizeResult``: output minimization routine.
    Get the result accessing `.x` attribute
    """

    def target_min_func(var_par: float):
        f_grid = func(var_par, *func_args)
        return energy_func_tht(f=f_grid, **energy_kwargs)

    return minimize_scalar(target_min_func, bounds=bounds, method="bounded")


def soliton_prototype(a: float, tht: np.ndarray) -> np.ndarray:
    """Generate a soliton in spherical polar grid points

    The soliton model is (1 - sin(tht)^a) e^(i * tanh(a*(tht - pi/2)) * pi/2)
    where it is narrower as larger is the parameter `a`. Due to the phase
    that must vary in pi from the south to north poles, the hyper tangent
    is used to smooth the variation.

    Return:
    -------
    ``numpy.ndarray``: soliton calculated in the grid points `tht`
    """
    sqrt_pi = np.sqrt(np.pi)
    norm_unit = np.sqrt(
        2
        + sqrt_pi * gamma(1 + a) / gamma(3 / 2 + a)
        - 2 * sqrt_pi * gamma(1 + a / 2) / gamma((3 + a) / 2)
    )
    return (
        (1 - np.sin(tht) ** a)
        * np.exp(1.0j * np.tanh((tht - np.pi / 2) * a) * np.pi / 2)
        / norm_unit
    )
