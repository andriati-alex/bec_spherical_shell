import os
import argparse
import bdg_driver
import numpy as np
from math import sqrt, pi


def run_from_imagtime(files_path, fname_prefix, m_vals, n_eigs):
    prefix = os.path.join(files_path, fname_prefix)
    suffix = "_imagtime.dat"
    eq_data = np.loadtxt(prefix + "_2species_equation" + suffix)
    if eq_data.ndim == 1:
        eq_data = eq_data.reshape(1, eq_data.size)
    njobs = eq_data.shape[0]
    tht_pts = int(eq_data[0, 0])
    frac_a = eq_data[0, 5]
    frac_b = eq_data[0, 6]
    vort_a = int(eq_data[0, 10])
    vort_b = int(eq_data[0, 11])
    ga_sweep = eq_data[:, 7] / (2 * pi)
    gb_sweep = eq_data[:, 8] / (2 * pi)
    gab_sweep = eq_data[:, 9] / (2 * pi)
    obs_data = np.loadtxt(prefix + "_2species_obs" + suffix)
    if obs_data.ndim == 1:
        obs_data = obs_data.reshape(1, obs_data.size)
    mu_a_sweep = obs_data[:, 2]
    mu_b_sweep = obs_data[:, 3]
    bdg = bdg_driver.BdGOperator(tht_pts, vort_a, vort_b, frac_a, frac_b)
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
        print("params {} {} {} {} {}".format(mu_a, mu_b, ga, gb, gab))
        for m in m_vals:
            eigs = bdg.lowlying_eig(m, mu_a, mu_b, ga, gb, gab, sa, sb)
            imag_eigs = np.sort(eigs.imag)[::-1]
            for i in range(n_eigs):
                out_file.write(" {:.5f}".format(imag_eigs[i]))
        out_file.write("\n")


def run_from_newton(files_path, fname_prefix, m_vals, n_eigs):
    prefix = os.path.join(files_path, fname_prefix)
    suffix = "_newton.dat"
    eq_data = np.loadtxt(prefix + "_equation" + suffix)
    if eq_data.ndim == 1:
        eq_data = eq_data.reshape(1, eq_data.size)
    obs_data = np.loadtxt(prefix + "_obs" + suffix)
    if obs_data.ndim == 1:
        obs_data = obs_data.reshape(1, obs_data.size)
    njobs = eq_data.shape[0]
    tht_pts = int(eq_data[0, 0])
    mu_a_sweep = eq_data[:, 3]
    mu_b_sweep = eq_data[:, 4]
    vort_a = int(eq_data[0, 5])
    vort_b = int(eq_data[0, 6])
    frac_a = obs_data[:, 5]
    frac_b = obs_data[:, 6]
    ga_sweep = obs_data[:, 7] / (2 * pi)
    gb_sweep = obs_data[:, 8] / (2 * pi)
    gab_sweep = obs_data[:, 9] / (2 * pi)
    norm_a = obs_data[:, 3]
    norm_b = obs_data[:, 4]
    params_pack = zip(
        mu_a_sweep,
        mu_b_sweep,
        ga_sweep,
        gb_sweep,
        gab_sweep,
        frac_a,
        frac_b,
        norm_a,
        norm_b,
    )
    out_file = open(prefix + "_stability.dat", "w")
    out_file.write("# {} largest imag part for each m\n".format(n_eigs))
    out_file.write("# m :")
    for m in m_vals:
        out_file.write(" {}".format(m))
    out_file.write("\n")
    for i, params_tuple in enumerate(params_pack):
        mu_a, mu_b, ga, gb, gab, fa, fb, na, nb = params_tuple
        bdg = bdg_driver.BdGOperator(tht_pts, vort_a, vort_b, fa, fb)
        job = i + 1
        print("Working on job [{:4d}/{}]".format(job, njobs))
        fname_a = prefix + "_speciesA" + "_job{}".format(job) + suffix
        fname_b = prefix + "_speciesB" + "_job{}".format(job) + suffix
        sa = np.loadtxt(fname_a) / sqrt(na)
        sb = np.loadtxt(fname_b) / sqrt(nb)
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
    methods = {"imagtime": run_from_imagtime, "newton": run_from_newton}
    p = argparse.ArgumentParser(
        usage="python %(prog)s `file_name_prefix` -mvals m1 ... mN ",
        description="Compute imag eigenvalues of several jobs",
    )
    p.add_argument(
        "fname_prefix",
        metavar="file_name_prefix",
        type=str,
        help="prefix for file names = string before `species` keyword",
    )
    p.add_argument(
        "-method",
        type=str,
        help="method used to compute the stationary states",
        default="imagtime",
        choices=methods.keys(),
    )
    p.add_argument(
        "--files-path",
        dest="files_path",
        type=str,
        help="path to .dat files",
        default=default_data_dir,
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
        default=3,
        type=int,
        help="number of eigenvalues to store for each case",
    )
    args = p.parse_args()
    kwargs = vars(args)
    func = methods[kwargs.pop("method")]
    func(**kwargs)
