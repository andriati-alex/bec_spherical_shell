import os
import argparse
import bdg_driver
import numpy as np
from math import pi, sqrt


def script_run(prefix, m_vals, n_eigs):
    suffix = "_imagtime.dat"
    eq_data = np.loadtxt(prefix + "_2species_equation" + suffix)
    njobs = eq_data.shape[0]
    tht_pts = int(eq_data[0, 0])
    frac_a = eq_data[0, 5]
    frac_b = eq_data[0, 6]
    vort_a = eq_data[0, 10]
    vort_b = eq_data[0, 11]
    ga_sweep = eq_data[:, 7] / (2 * pi)
    gb_sweep = eq_data[:, 8] / (2 * pi)
    gab_sweep = eq_data[:, 9] / (2 * pi)
    obs_data = np.loadtxt(prefix + "_2species_obs" + suffix)
    mu_a_sweep = obs_data[:, 2]
    mu_b_sweep = obs_data[:, 3]
    fallback_a = np.sqrt(0.5) * np.ones(tht_pts)
    fallback_b = np.sqrt(0.5) * np.ones(tht_pts)
    bdg = bdg_driver.BdGOperator(
        fallback_a, fallback_b, vort_a, vort_b, frac_a, frac_b
    )
    params_pack = zip(mu_a_sweep, mu_b_sweep, ga_sweep, gb_sweep, gab_sweep)
    out_file = open(prefix + "_stability.dat", "w")
    out_file.write("# {} largest imag part for each m\n".format(n_eigs))
    out_file.write("# m :\n")
    for m in m_vals:
        out_file.write(" {}".format(m))
    out_file.write("\n")
    for i, (mu_a, mu_b, ga, gb, gab) in enumerate(params_pack):
        job = i + 1
        print("Working on job [{}/{}]".format(job, njobs))
        fname_a = prefix + "_speciesA" + "_job{}".format(job) + suffix
        fname_b = prefix + "_speciesB" + "_job{}".format(job) + suffix
        raw_sa = np.loadtxt(fname_a, dtype=np.complex128)
        raw_sb = np.loadtxt(fname_b, dtype=np.complex128)
        sa = (
            raw_sa * np.exp(-1.0j * np.arctan2(raw_sa.imag, raw_sa.real))
        ).real
        sb = (
            raw_sb * np.exp(-1.0j * np.arctan2(raw_sb.imag, raw_sb.real))
        ).real
        for m in m_vals:
            eigs = bdg.lowlying_eig(m, mu_a, mu_b, ga, gb, gab, sa, sb)
            imag_eigs = np.sort(eigs.imag)[::-1]
            for i in range(n_eigs):
                out_file.write(" {:.5f}".format(imag_eigs[i]))
        out_file.write("\n")


if __name__ == "__main__":
    default_data_dir = os.path.join(
        os.path.expanduser("~"), "programs/bec_spherical_shell/output"
    )
    p = argparse.ArgumentParser(
        usage="python %(prog)s `files_prefix` ",
        description="Compute imag eigenvalues of several jobs",
    )
    p.add_argument(
        "prefix",
        metavar="files_prefix",
        type=str,
        help="prefix for file names = string before `species` keyword",
    )
    p.add_argument(
        "-mvals",
        dest="m_vals",
        nargs="+",
        type=int,
        help="azimuthal number of pertubation",
    )
    p.add_argument(
        "-neigs",
        dest="n_eigs",
        default=5,
        type=int,
        help="number of eigenvalues to store for each case",
    )
    args = p.parse_args()
    args.prefix = os.path.join(default_data_dir, args.prefix)
    script_run(**vars(args))
