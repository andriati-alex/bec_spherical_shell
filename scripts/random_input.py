import os
import argparse
from random import randint
import numpy as np
from math import pi, sqrt, factorial
from scipy.integrate import simps
from scipy.special import sph_harm


def general_spherical_state(lmax, phi_pts, tht_pts, seed):
    if not isinstance(seed, int):
        raise ValueError("Seed provided {} is not integer".format(seed))
    rand_generator = np.random.default_rng(seed=seed)
    dphi = 2 * pi / (phi_pts - 1)
    dtht = pi / (tht_pts - 1)
    phi, tht = np.meshgrid(
        np.linspace(0, 2 * pi, phi_pts), np.linspace(0, pi, tht_pts)
    )
    grid_noise = np.zeros([tht_pts, phi_pts], dtype=np.complex128)
    grid_noise[1 : tht_pts - 1, 1 : phi_pts - 1] = (
        0.25
        * (
            rand_generator.random([tht_pts - 2, phi_pts - 2])
            + 1.0j * rand_generator.random([tht_pts - 2, phi_pts - 2])
            - 0.5
            - 0.5j
        )
        * np.sin(tht)
    )
    psi = sph_harm(0, 0, phi, tht) + grid_noise
    number_sph_harm = int((2 * lmax + 1 + 1) * (lmax + 1) / 2)
    w = 5 * (
        rand_generator.random(number_sph_harm)
        + 1.0j * rand_generator.random(number_sph_harm)
        - 0.5
        - 0.5j
    )
    k = 0
    for l in range(1, lmax + 1):
        for m in range(-l, l + 1):
            psi = psi + w[k] * sph_harm(m, l, phi, tht) / sqrt(factorial(l))
            k = k + 1
    abs_square = abs(psi) ** 2
    integral_phi = np.zeros(psi.shape[0], dtype=np.float64)
    for i in range(psi.shape[0]):
        integral_phi[i] = simps(abs_square[i], dx=dphi)
    return psi / sqrt(simps(np.sin(tht[:, 0]) * integral_phi, dx=dtht))


if __name__ == "__main__":
    default_input_dir = os.path.join(
        os.path.expanduser("~"), "programs/bec_spherical_shell/input"
    )
    p = argparse.ArgumentParser(
        usage="python %(prog)s n_phi n_theta dt n_dt [optional_args] ",
        description="Initial condition generator of random"
        " quantum states in a spherical shell",
    )
    p.add_argument(
        "phi_pts",
        metavar="n_phi",
        action="store",
        type=int,
        help="number of grid points in `phi` axis",
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
        default=15,
        type=int,
        help="maximum spherical harmonic number l",
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

    psi_a = general_spherical_state(
        args.lmax, args.phi_pts, args.tht_pts, args.seed_a
    )
    psi_b = general_spherical_state(
        args.lmax, args.phi_pts, args.tht_pts, args.seed_b
    )

    # Record initial state A
    np.savetxt(
        os.path.join(args.input_dir, "rnd_state_speciesA_init.dat"),
        psi_a.reshape(args.phi_pts * args.tht_pts).T,
        fmt="%.15E",
    )

    # Record initial state B
    np.savetxt(
        os.path.join(args.input_dir, "rnd_state_speciesB_init.dat"),
        psi_b.reshape(args.phi_pts * args.tht_pts).T,
        fmt="%.15E",
    )

    # Configure domain
    f = open(os.path.join(args.input_dir, "rnd_state_domain.dat"), "w")
    f.write(
        "{} {} {} {}\n".format(args.phi_pts, args.tht_pts, args.dt, args.ndt)
    )
    f.close()

    # Default values are set for equation parameters
    f = open(os.path.join(args.input_dir, "rnd_state_eq.dat"), "w")
    f.write(
        "-0.5 0.0 {:.2f} {:.2f} {:.5f} {:.5f} {:.5f}\n".format(
            args.frac_a, args.frac_b, args.ga, args.gb, args.gab
        )
    )
    f.close()
