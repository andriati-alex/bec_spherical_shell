import os
import argparse
from random import randint
import numpy as np
from math import pi, sqrt, factorial
from scipy.integrate import simps


def theta_random_state(lmax, tht_pts, azi, seed):
    if not isinstance(seed, int):
        raise ValueError("Seed provided {} is not integer".format(seed))
    if not isinstance(azi, int):
        raise ValueError("azimuthal number {} is not integer".format(azi))
    rand_generator = np.random.default_rng(seed=seed)
    dtht = pi / (tht_pts - 1)
    tht = np.linspace(0, pi, tht_pts)
    grid_noise = np.zeros(tht_pts, dtype=np.complex128)
    grid_noise[1 : tht_pts - 1] = 0.1 * (
        rand_generator.random([tht_pts - 2])
        + 1.0j * rand_generator.random([tht_pts - 2])
        - 0.5
        - 0.5j
    )
    if azi == 0:
        trig_func = np.cos
        psi = np.ones(tht_pts, dtype=np.complex128) + grid_noise
    else:
        trig_func = np.sin
        psi = np.sin(tht) + grid_noise
    w = 5 * (
        rand_generator.random(lmax)
        + 1.0j * rand_generator.random(lmax)
        - 0.5
        - 0.5j
    )
    for j in range(1, lmax + 1):
        psi = psi + w[j - 1] * trig_func(j * tht) / sqrt(factorial(j))
    return psi / sqrt(simps(np.sin(tht) * abs(psi) ** 2, dx=dtht))


if __name__ == "__main__":
    default_input_dir = os.path.join(
        os.path.expanduser("~"), "projects/bec_spherical_shell/input"
    )
    p = argparse.ArgumentParser(
        usage="python %(prog)s n_theta dt n_dt [optional_args] ",
        description="Initial condition generator of random quantum"
        " states in a spherical shell with fixed azimuthal number",
    )
    p.add_argument(
        "tht_pts",
        metavar="n_theta",
        action="store",
        type=int,
        help="number of grid points in `theta` axis",
    )
    p.add_argument(
        "dt",
        action="store",
        type=float,
        help="time step size",
    )
    p.add_argument(
        "ndt",
        metavar="n_dt",
        action="store",
        type=int,
        help="number of time steps to propagate initial condition",
    )
    p.add_argument(
        "-path",
        dest="input_dir",
        action="store",
        default=default_input_dir,
        type=str,
        help="path to save files generated",
    )
    p.add_argument(
        "-lmax",
        action="store",
        default=45,
        type=int,
        help="number of states in basis expansion with random coefficients",
    )
    p.add_argument(
        "--seed-a",
        action="store",
        dest="seed_a",
        default=randint(1, 1000000),
        type=int,
        help="Integer seed to generate random numbers",
    )
    p.add_argument(
        "--seed-b",
        action="store",
        dest="seed_b",
        default=randint(1, 1000000),
        type=int,
        help="Integer seed to generate random numbers",
    )
    p.add_argument(
        "--azi-a",
        action="store",
        dest="azi_a",
        default=0,
        type=int,
        help="azimuthal number of phi exponential",
    )
    p.add_argument(
        "--azi-b",
        action="store",
        dest="azi_b",
        default=0,
        type=int,
        help="azimuthal number of phi exponential",
    )
    p.add_argument(
        "--frac-a",
        action="store",
        dest="frac_a",
        default=0.5,
        type=float,
        help="fraction of atoms of species A",
    )
    p.add_argument(
        "--frac-b",
        action="store",
        dest="frac_b",
        default=0.5,
        type=float,
        help="fraction of atoms of spceies B",
    )
    p.add_argument(
        "-ga",
        action="store",
        dest="ga",
        default=10,
        type=float,
        help="interaction strength species A",
    )
    p.add_argument(
        "-gb",
        action="store",
        dest="gb",
        default=10,
        type=float,
        help="interaction strength species B",
    )
    p.add_argument(
        "-gab",
        action="store",
        dest="gab",
        default=5,
        type=float,
        help="interspecies interaction strength",
    )
    args = p.parse_args()

    if args.input_dir == default_input_dir:
        os.makedirs(args.input_dir, exist_ok=True)
    else:
        if not os.path.isdir(args.input_dir):
            dir_msg = "path {} does not exist, create it [y/n]? ".format(
                args.input_dir
            )
            must_create = input(dir_msg) == "y"
            if must_create:
                os.makedirs(args.input_dir)
            else:
                print("process aborted")
                exit()

    psi_a = theta_random_state(
        args.lmax, args.tht_pts, args.azi_a, args.seed_a
    )
    psi_b = theta_random_state(
        args.lmax, args.tht_pts, args.azi_b, args.seed_b
    )

    # Record initial state A
    np.savetxt(
        os.path.join(args.input_dir, "rndtheta_state_speciesA_init.dat"),
        psi_a.T,
        fmt="%.15E",
    )

    # Record initial state B
    np.savetxt(
        os.path.join(args.input_dir, "rndtheta_state_speciesB_init.dat"),
        psi_b.T,
        fmt="%.15E",
    )

    # Configure domain
    f = open(os.path.join(args.input_dir, "rndtheta_state_domain.dat"), "w")
    f.write("{} {} {}\n".format(args.tht_pts, args.dt, args.ndt))
    f.close()

    # Default values are set for equation parameters
    f = open(os.path.join(args.input_dir, "rndtheta_state_eq.dat"), "w")
    f.write(
        "-0.5 0.0 {:.2f} {:.2f} {:.5f} {:.5f} {:.5f} {:d} {:d}\n".format(
            args.frac_a,
            args.frac_b,
            args.ga,
            args.gb,
            args.gab,
            args.azi_a,
            args.azi_b,
        )
    )
    f.close()
