import os
import argparse
import random
import numpy as np
from math import pi, sqrt, factorial
from scipy.integrate import simps
from scipy.special import sph_harm
from constrained_minimization import StateSphericalBasis


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
        "-lmin",
        action="store",
        default=0,
        type=int,
        help="maximum spherical harmonic number l",
    )
    p.add_argument(
        "-lmax",
        action="store",
        default=4,
        type=int,
        help="maximum spherical harmonic number l",
    )
    p.add_argument(
        "--m-min",
        action="store",
        dest="m_min",
        default=0,
        type=int,
        help="maximum spherical harmonic number l",
    )
    p.add_argument(
        "-seed",
        action="store",
        default=random.randint(1, 100000),
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

    state_in_basis = StateSphericalBasis(
        args.lmax, args.lmin, args.m_min, args.phi_pts, args.tht_pts, args.seed
    )
    e, psi_a, psi_b = state_in_basis.minimize_energy(
        args.frac_a, args.frac_b, args.ga, args.gb, args.gab
    )
    print("Initial condition energy : {:.3f}".format(e))

    # Record initial state A
    np.savetxt(
        os.path.join(args.input_dir, "opt_state_speciesA_init.dat"),
        psi_a.reshape(args.phi_pts * args.tht_pts).T,
        fmt="%.15E",
    )

    # Record initial state B
    np.savetxt(
        os.path.join(args.input_dir, "opt_state_speciesB_init.dat"),
        psi_b.reshape(args.phi_pts * args.tht_pts).T,
        fmt="%.15E",
    )

    # Configure domain
    f = open(os.path.join(args.input_dir, "opt_state_domain.dat"), "w")
    f.write(
        "{} {} {} {}\n".format(args.phi_pts, args.tht_pts, args.dt, args.ndt)
    )
    f.close()

    # Default values are set for equation parameters
    f = open(os.path.join(args.input_dir, "opt_state_eq.dat"), "w")
    f.write(
        "-0.5 0.0 {:.2f} {:.2f} {:.5f} {:.5f} {:.5f}\n".format(
            args.frac_a, args.frac_b, args.ga, args.gb, args.gab
        )
    )
    f.close()
