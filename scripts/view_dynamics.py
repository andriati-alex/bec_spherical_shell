import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import pi
from mayavi import mlab


def export_phase_colorbar(fname="phase_colorbar.pdf", fig_height=1.75):
    mpl.rcParams["font.size"] = 9
    mpl.rcParams["text.usetex"] = True
    fig = plt.figure(figsize=(0.4, fig_height))
    ax = fig.add_axes([0.4, 0.1, 0.2, 0.8])
    cmap_norm = mpl.colors.Normalize(vmin=0, vmax=2 * pi)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(cmap=mpl.cm.hsv, norm=cmap_norm), cax=ax
    )
    ax.tick_params(axis="both", pad=0)
    cbar.set_ticks([0, pi, 2 * pi])
    cbar.set_ticklabels(["$0$", "$\pi$", "$2\pi$"])
    cbar.ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
    plt.savefig(fname, dpi=512, bbox_inches="tight")


def export_density_colorbar(
    fname="density_colorbar.pdf", den_min=0, den_max=0, fig_height=1.75
):
    mpl.rcParams["font.size"] = 9
    mpl.rcParams["text.usetex"] = True
    fig = plt.figure(figsize=(0.4, fig_height))
    ax = fig.add_axes([0.4, 0.1, 0.2, 0.8])
    cmap_norm = mpl.colors.Normalize(vmin=den_min, vmax=den_max)
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=mpl.cm.jet, norm=cmap_norm), cax=ax
    )
    ax.tick_params(axis="both", pad=0)
    plt.savefig(fname, dpi=512, bbox_inches="tight")


class TimeDependentState:
    """
    Class to handle time dependent states from .dat files
    Most part of functionality is devoted to get density
    arrays for a given time frame index. See the function
    `frame_from_time` to obtain the index of the closest
    frame of a certain time to help in the input of other
    functions

    Main attributes
    ---
    `self.time_den` : ``numpy.array``
        time instants the states were recorded
    `self.time` : ``numpy.array``
        all time instants of numerical time evolution
        also time instants the observables were recorded
    `self.overlap` : ``numpy.array``
        overlap functional for all time instants
    `self.lz(a/b)` : ``numpy.array``
        angular momentum for all time instants
    """

    def __init__(self, data_path, prefix):
        """
        Parameters
        ---
        `data_path` : ``str``
            valid path to a folder with all .dat files from main program output
        `prefix` : ``str``
            file prefix convention

        """
        file_prefix = os.path.join(data_path, prefix)
        eq_filename = file_prefix + "_2species_equation_realtime.dat"
        if not os.path.isfile(eq_filename):
            raise IOError("\nEquation file {} not found\n".format(eq_filename))
        eq_data = np.loadtxt(file_prefix + "_2species_equation_realtime.dat")
        self.phi = np.linspace(0, 2 * pi, int(eq_data[0]))
        self.theta = np.linspace(0, pi, int(eq_data[1]))
        self.frac_a = eq_data[6]
        self.frac_b = eq_data[7]
        self.ntheta = self.theta.size
        self.nphi = self.phi.size
        self.a_fname = file_prefix + "_speciesA_realtime.dat"
        self.b_fname = file_prefix + "_speciesB_realtime.dat"
        self.time_den = np.loadtxt(file_prefix + "_log_realtime.dat")[:, 0]
        self.den_max = self.__all_time_den_max()
        self.den_min = self.__all_time_den_min()
        obs_data = np.loadtxt(file_prefix + "_obs_realtime.dat")
        self.time = obs_data[:, 0]
        self.overlap = obs_data[:, 1]
        self.lza = obs_data[:, 2]
        self.lzb = obs_data[:, 3]

    def frame_from_time(self, t):
        """ Get the closest frame index of a given time instant `t` """
        return abs(self.time_den - t).argmin()

    def frame_min_overlap(self):
        """ Get frame with minimum overlap along time evolution """
        return self.frame_from_time(self.time[self.overlap.argmin()])

    def __all_time_den_min(self):
        dena = abs(np.loadtxt(self.a_fname)) ** 2
        denb = abs(np.loadtxt(self.b_fname)) ** 2
        den = self.frac_a * dena + self.frac_b * denb
        return min(dena, denb, den)

    def __all_time_den_max(self):
        dena = abs(np.loadtxt(self.a_fname)) ** 2
        denb = abs(np.loadtxt(self.b_fname)) ** 2
        den = self.frac_a * dena + self.frac_b * denb
        return max(dena, denb, den)

    def get_states(self, frame_index=0):
        """ Get ``tuple`` with the two states in a given `frame_index` """
        sa = np.loadtxt(self.a_fname, skip_rows=frame_index).reshape(
            self.ntheta, self.nphi
        )
        sb = np.loadtxt(self.b_fname, skip_rows=frame_index).reshape(
            self.ntheta, self.nphi
        )
        return sa, sb

    def get_densities(self, frame_index=0):
        """ Return ``tuple`` with the three densities in the `frame_index` """
        sa, sb = self.get_states(frame_index)
        dena = abs(sa) ** 2
        denb = abs(sb) ** 2
        return dena, denb, self.frac_a * dena + self.frac_b * denb

    def cartesian_grid(self, radius=1, offset=(0, 0, 0)):
        """
        Generate cartesian grid points of sphere centered at `offset`

        Parameters
        ---
        `radius` : ``float``
            sphere radius
        `offset` :  ``tuple`` (float, float, float)
            cartesian coordinates of sphere center

        """
        phi_grid, tht_grid = np.meshgrid(self.phi, self.theta)
        x = radius * np.sin(tht_grid) * np.cos(phi_grid) + offset[0]
        y = radius * np.sin(tht_grid) * np.sin(phi_grid) + offset[1]
        z = radius * np.cos(tht_grid) + offset[2]
        return (x, y, z)

    def normalize_density(self, den):
        """
        Normalize density in [0, 1] interval considering
        the all time maximum values

        """
        return (den - self.den_min) / (self.den_max - self.den_min)

    def max_density_angles(self, den):
        """
        Return (``float``, ``float``) with phi, theta correpondingo to
        maximum density

        """
        stride = den.argmax()
        i, j = stride // self.theta.size, stride % self.phi.size
        return self.phi[j], self.theta[i]

    def frame_density(self, frame_index=0, view_angle=(0, 0), cmap="jet"):
        """
        Plot 3 spheres with colormaps corresponding to densities
        From left to right : density A / density B / overall density

        Parmeters
        ---
        `frame_index` : ``int``
            number of recorded frame in .dat files
        `view_angle` : ``tuple`` (float, float)
            angles of mayavi view. See documentation

        """
        mlab.figure(
            bgcolor=(1.0, 1.0, 1.0),
            size=(1100, 600),
        )
        mlab.clf()
        xa, ya, za = self.cartesian_grid(offset=(-2.15, 0, 0))
        xb, yb, zb = self.cartesian_grid(offset=(0, 0, 0))
        x, y, z = self.cartesian_grid(offset=(2.15, 0, 0))
        den_a, den_b, den = self.get_densities(self.frame_min_overlap())
        norm_den_a = self.normalize_density(den_a)
        norm_den_b = self.normalize_density(den_b)
        norm_total_den = self.frac_a * norm_den_a + self.frac_b * norm_den_b
        mlab.mesh(xa, ya, za, scalars=norm_den_a, colormap=cmap)
        mlab.mesh(xb, yb, zb, scalars=norm_den_b, colormap=cmap)
        mlab.mesh(x, y, z, scalars=norm_total_den, colormap=cmap)
        mlab.view(*view_angle)
