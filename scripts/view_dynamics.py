import os
import linecache
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from math import pi
from mayavi import mlab


def gen_specific_row(filepath, row_index=0):
    """ yield linecache.getline(filepath, row_index + 1) """
    yield linecache.getline(filepath, row_index + 1)


def export_phase_colorbar(fname="phase_cbar.pdf", fig_height=1.75):
    """ Save pdf file with colorbar ranging from 0 to 2 * pi """
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
    fname="density_cbar.pdf", vmin=0, vmax=1, fig_height=1.75
):
    """ Save pdf file with colorbar given min and max values """
    mpl.rcParams["font.size"] = 9
    mpl.rcParams["text.usetex"] = True
    fig = plt.figure(figsize=(0.4, fig_height))
    ax = fig.add_axes([0.4, 0.1, 0.2, 0.8])
    cmap_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
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
        time instants the full wave functions were recorded
    `self.time_obs` : ``numpy.array``
        all time instants of numerical time evolution
        time instants some observables were recorded
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
            path to a folder with all .dat files from main program output
        `prefix` : ``str``
            file prefix convention (common part for all .dat file names)

        """
        file_prefix = os.path.join(data_path, prefix)
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
        obs_data = np.loadtxt(file_prefix + "_obs_realtime.dat")
        self.time_obs = obs_data[:, 0]
        self.overlap = obs_data[:, 1]
        self.lza = obs_data[:, 2]
        self.lzb = obs_data[:, 3]

    def frame_from_time(self, t):
        """ Get the closest frame index of a given time instant `t` """
        return abs(self.time_den - t).argmin()

    def frame_min_overlap(self):
        """ Get frame index with minimum overlap along time evolution """
        return self.frame_from_time(self.time_obs[self.overlap.argmin()])

    def nearest_theta_index(self, theta_value):
        """ Return index of theta array with value closest to `theta_value` """
        return abs(self.theta - theta_value).argmin()

    def get_state_a(self, frame_index=0):
        """ Return state of species A at some `frame_index` """
        return np.loadtxt(
            gen_specific_row(self.a_fname, frame_index), dtype=np.complex128
        ).reshape(self.ntheta, self.nphi)

    def get_state_b(self, frame_index=0):
        """ Return state of species B at some `frame_index` """
        return np.loadtxt(
            gen_specific_row(self.b_fname, frame_index), dtype=np.complex128
        ).reshape(self.ntheta, self.nphi)

    def get_states(self, frame_index=0):
        """ Get ``tuple`` with the two states at given `frame_index` """
        return self.get_state_a(frame_index), self.get_state_b(frame_index)

    def get_density_a(self, frame_index=0):
        """ Return density of species A at given `frame_index` """
        return abs(self.get_state_a(frame_index)) ** 2

    def get_density_b(self, frame_index=0):
        """ Return density of species B at given `frame_index` """
        return abs(self.get_state_b(frame_index)) ** 2

    def get_densities(self, frame_index=0):
        """ Return ``tuple`` with A,B densities at `frame_index` """
        return self.get_density_a(frame_index), self.get_density_b(frame_index)

    def cartesian_grid(self, radius=1, offset=(0, 0, 0)):
        """
        Generate cartesian grid points of surface in
        spherical coordinates centered at `offset`

        Parameters
        ---
        `radius` : ``float`` or ``numpy.array``
            surface radius as function of spherical angles
        `offset` :  ``tuple`` (float, float, float)
            cartesian coordinates of surface center (radius origin)

        Return
        ---
        ``tuple``
            (`x`, `y`, `z`) numpy grids with surface surface coordinates

        """
        phi_grid, tht_grid = np.meshgrid(self.phi, self.theta)
        x = radius * np.sin(tht_grid) * np.cos(phi_grid) + offset[0]
        y = radius * np.sin(tht_grid) * np.sin(phi_grid) + offset[1]
        z = radius * np.cos(tht_grid) + offset[2]
        return (x, y, z)

    def normalize_density(self, den):
        """ Normalize density in [0, 1] interval """
        return (den - den.min()) / (den.max() - den.min())

    def max_density_angles(self, den):
        """
        Return tuple of spherical angles (phi, theta)
        corresponding to point of maximum density
        """
        stride = den.argmax()
        i, j = stride // self.theta.size, stride % self.phi.size
        return self.phi[j], self.theta[i]

    def density_lines(self, theta_ind, frame_index=0):
        """
        Return tuple of 2 arrays with density along phi for fixed theta

        Parameters
        ---
        `theta_ind` : ``int``
            index of theta array values (use self.nearest_theta_index to help)
        `frame_index` : ``int``
            index of the frame in time evolution

        Return
        ---
        ``(numpy.array, numpy.array)``
            density species A, density of species B

        """
        dena, denb = self.get_densities(frame_index)
        return dena[theta_ind], denb[theta_ind]

    def plot_all_densities(self, frame_index=0, view_args=(0, 0), cmap="jet"):
        """
        Plot 3 spheres with colormaps corresponding to densities
        From left to right : density A / density B / overall density

        Parmeters
        ---
        `frame_index` : ``int``
            frame number (row in .dat files)
        `view_args` : ``tuple`` (float, float)
            mayavi view arguments. See documentation
            https://docs.enthought.com/mayavi/mayavi
        `cmap` : ``str``
            available colormap name. Consult matplotlib API

        """
        mlab.figure(bgcolor=(1.0, 1.0, 1.0), size=(1100, 600))
        mlab.clf()
        xa, ya, za = self.cartesian_grid(offset=(-2.15, 0, 0))
        xb, yb, zb = self.cartesian_grid(offset=(0, 0, 0))
        x, y, z = self.cartesian_grid(offset=(2.15, 0, 0))
        den_a, den_b = self.get_densities(frame_index)
        norm_den_a = self.normalize_density(den_a)
        norm_den_b = self.normalize_density(den_b)
        norm_total_den = self.frac_a * norm_den_a + self.frac_b * norm_den_b
        mlab.mesh(xa, ya, za, scalars=norm_den_a, colormap=cmap)
        mlab.mesh(xb, yb, zb, scalars=norm_den_b, colormap=cmap)
        mlab.mesh(x, y, z, scalars=norm_total_den, colormap=cmap)
        mlab.view(*view_args)

    def plot_density(
        self, frame_index=0, species="A", view_args=(0, 0), cmap="jet"
    ):
        """
        Plot and show density mapped to colors in a sphere

        Parmeters
        ---
        `frame_index` : ``int``
            frame number (row in .dat files)
        `species` : "A" or "B"
            species type according to filename convention
        `view_args` : ``tuple`` (float, float)
            mayavi view arguments. See documentation
            https://docs.enthought.com/mayavi/mayavi
        `cmap` : ``str``
            available colormap name. Consult matplotlib API

        """
        mlab.figure(bgcolor=(1.0, 1.0, 1.0), size=(800, 600))
        mlab.clf()
        x, y, z = self.cartesian_grid()
        if species.upper() == "A":
            norm_den = self.normalize_density(self.get_density_a(frame_index))
        else:
            norm_den = self.normalize_density(self.get_density_b(frame_index))
        mlab.mesh(x, y, z, scalars=norm_den, colormap=cmap)
        mlab.view(*view_args)

    def plot_radial_surface(
        self, frame_index=0, species="A", view_args=(0, 0)
    ):
        """
        Use density as radial distance and phase as colors

        Parmeters
        ---
        `frame_index` : ``int``
            frame number (row in .dat files)
        `species` : "A" or "B"
            species type according to filename convention
        `view_args` : ``tuple`` (float, float)
            mayavi view arguments. See documentation
            https://docs.enthought.com/mayavi/mayavi

        """
        mlab.figure(bgcolor=(1.0, 1.0, 1.0), size=(800, 600))
        mlab.clf()
        if species.upper() == "A":
            state = self.get_state_a(frame_index)
            norm_den = self.normalize_density(self.get_density_a(frame_index))
        else:
            state = self.get_state_b(frame_index)
            norm_den = self.normalize_density(self.get_density_b(frame_index))
        phase = np.arctan2(state.imag, state.real)
        x, y, z = self.cartesian_grid(norm_den)
        mlab.mesh(x, y, z, scalars=phase, colormap="hsv")
        mlab.view(*view_args)

    def __videoclip_frame_radial_surface(self, t):
        """
        Return raw image data in current animation time `anim_t`
        This is an auxiliar method used in `self.videoclip_obj`
        """
        mlab.clf()
        prog_ratio = t / self.__anim_duration
        data_period = self.__t_stop - self.__t_start
        t_scale = data_period * prog_ratio + self.__t_start
        frame_index = self.frame_from_time(t_scale)
        state = self.get_state_a(frame_index)
        phase = np.arctan2(state.imag, state.real)
        den = self.normalize_density(self.get_density_a(frame_index))
        x, y, z = self.cartesian_grid(den)
        mlab.mesh(x, y, z, scalars=phase, colormap="hsv")
        mlab.view(*(self.__anim_view_args))
        return mlab.screenshot(antialiased=True)

    def __videoclip_frame(self, t):
        """
        Return raw image data in current animation time `anim_t`
        This is an auxiliar method used in `self.videoclip_obj`
        """
        mlab.clf()
        prog_ratio = t / self.__anim_duration
        data_period = self.__t_stop - self.__t_start
        t_scale = data_period * prog_ratio + self.__t_start
        frame_index = self.frame_from_time(t_scale)
        den = self.normalize_density(self.get_density_a(frame_index))
        mlab.mesh(*(self.__anim_sphere_grid), scalars=den, colormap="jet")
        mlab.view(*(self.__anim_view_args))
        return mlab.screenshot(antialiased=True)

    def videoclip_obj(
        self,
        start=0,
        stop=None,
        anim_duration=30,
        fig_size=(800, 600),
        view_args=(),
        display_phase=False,
        end_rest=0.5,
    ):
        """Generate animation object using `moviepy` module
        Consult methods to write the returned videoclip

        Parameters
        ---
        `start` : ``float``
            initial time in system units. Must be in self.time_den
        `stop` : ``float``
            Final time in system units. Must be in self.time_den
        `anim_duration` : ``int``
            output movie duration in seconds
        `fig_size` : ``tuple`` (width, height)
            figure size to display the animation
        `view_args` : ``tuple`` of floats
            mayavi view arguments. See documentation
            https://docs.enthought.com/mayavi/mayavi
        `display_phase` : ``bool``
            Whether to use radial plot and display phase as colormap (True)
        `end_rest` : ``float``
            Time in seconds to freeze in the last frame

        Return
        ---
        ``moviepy.editor.VideoClip``
            video clip animation object

        """
        mlab.figure(bgcolor=(1.0, 1.0, 1.0), size=fig_size)
        mlab.clf()
        mlab.view(*view_args)
        f = mlab.gcf()
        f.scene._lift()
        if stop is None:
            stop = self.time_den[-1]
        self.__t_start = start
        self.__t_stop = stop
        self.__anim_view_args = view_args
        self.__anim_sphere_grid = self.cartesian_grid()
        self.__anim_duration = anim_duration - end_rest
        if display_phase:
            anim_method = self.__videoclip_frame_radial_surface
        else:
            anim_method = self.__videoclip_frame
        return mpy.VideoClip(anim_method, duration=anim_duration)
