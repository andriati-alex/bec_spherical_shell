import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
from numpy import pi, sqrt

mpl.rcParams["text.usetex"] = True
mpl.rcParams["figure.subplot.hspace"] = 0.0
mpl.rcParams["figure.subplot.wspace"] = 0.0


class DensityPlot:
    def __init__(self, data_path, prefix, view_tool="matplotlib"):
        file_prefix = os.path.join(data_path, prefix)
        eq_filename = file_prefix + "_2species_equation_imagtime.dat"
        if not os.path.isfile(eq_filename):
            raise IOError("\nEquation file {} not found\n".format(eq_filename))
        eq_data = np.loadtxt(file_prefix + "_2species_equation_imagtime.dat")
        self.phi = np.linspace(0, 2 * np.pi, int(eq_data[0]))
        self.theta = np.linspace(0, np.pi, int(eq_data[1]))
        self.nabla_coef = eq_data[4]
        self.rotation_freq = eq_data[5]
        self.frac_a = eq_data[6]
        self.frac_b = eq_data[7]
        self.inter_a = eq_data[8]
        self.inter_b = eq_data[9]
        self.inter_ab = eq_data[10]
        self.state_a = np.loadtxt(
            file_prefix + "_speciesA_imagtime.dat", dtype=np.complex128
        ).reshape(self.theta.size, self.phi.size)
        self.state_b = np.loadtxt(
            file_prefix + "_speciesB_imagtime.dat", dtype=np.complex128
        ).reshape(self.theta.size, self.phi.size)
        self.abs_square_a = abs(self.state_a) ** 2
        self.abs_square_b = abs(self.state_b) ** 2
        self.view_tool = view_tool

    def show_in_one_sphere(self):
        fig = plt.figure(figsize=plt.figaspect(1.0))
        ax = fig.gca(projection="3d")
        colors_a = self.__rgba_colors("red", self.abs_square_a)
        colors_b = self.__rgba_colors("blue", self.abs_square_b)
        surf_a = self.__plot_rgba_sphere(ax, 1.0, colors_a)
        surf_b = self.__plot_rgba_sphere(ax, 1.0, colors_b)
        ax.set_axis_off()
        self.__set_view_angle_mpl(ax)
        plt.show()

    def show_in_two_spheres(self):
        # create figure
        fig = plt.figure(figsize=(8, 6))
        axa = fig.add_subplot(121, projection="3d")
        axb = fig.add_subplot(122, projection="3d")
        self.__plot_cmap_sphere(fig, axa, self.abs_square_a)
        self.__plot_cmap_sphere(fig, axb, self.abs_square_b)
        axa.set_title("$|\Psi_1(\\theta,\phi)|^2$")
        axb.set_title("$|\Psi_2(\\theta,\phi)|^2$")
        plt.show()
        # plt.savefig("teste.png", dpi=1024, bbox_inches="tight")

    def show_mayavi_oriented(self):
        mlab.figure(
            1,
            bgcolor=(1.0, 1.0, 1.0),
            fgcolor=(0.0, 0.0, 0.0),
            size=(1100, 600),
        )
        mlab.clf()
        phi_den_max, theta_den_max = self.__max_density_angles()
        max_den_direction = np.array(
            [
                np.sin(theta_den_max) * np.cos(phi_den_max),
                np.sin(theta_den_max) * np.sin(phi_den_max),
                np.cos(theta_den_max),
            ]
        )
        self.__set_view_angle_mayavi()
        self.__mayavi_sphere(self.abs_square_a)
        camera_position, focus_position = mlab.move()
        cross_vec = np.cross(camera_position, max_den_direction)
        offset_second_sphere = (
            2.5 * cross_vec / sqrt((abs(cross_vec) ** 2).sum())
        )
        self.__mayavi_sphere(self.abs_square_b, offset=offset_second_sphere)
        mlab.show()

    def show_mayavi_raw(self):
        mlab.figure(
            1,
            bgcolor=(1.0, 1.0, 1.0),
            fgcolor=(0.0, 0.0, 0.0),
            size=(1100, 600),
        )
        mlab.clf()
        mlab.view(azimuth=pi, elevation=pi / 2)
        self.__mayavi_sphere(self.abs_square_a)
        self.__mayavi_sphere(self.abs_square_b, offset=(2.5, 0.0, 0.0))
        mlab.show()

    def __mayavi_sphere(self, den, radius=1, offset=(0, 0, 0), cmap="jet"):
        x, y, z = self.__cartesian_grid(radius, offset)
        mlab.mesh(x, y, z, scalars=den, colormap=cmap)

    def __cartesian_grid(self, radius=1, offset=(0, 0, 0)):
        phi_grid, tht_grid = np.meshgrid(self.phi, self.theta)
        # sphere of unit radius
        x = radius * np.sin(tht_grid) * np.cos(phi_grid) + offset[0]
        y = radius * np.sin(tht_grid) * np.sin(phi_grid) + offset[1]
        z = radius * np.cos(tht_grid) + offset[2]
        return (x, y, z)

    def __density_contrast(self, den):
        return 1.0 - den.min() / den.max()

    def __normalize_density(self, den):
        return (den - den.min()) / (den.max() - den.min())

    def __max_density_angles(self):
        stride = self.abs_square_a.argmax()
        i, j = int(stride / self.theta.size), stride % self.phi.size
        return self.phi[j], self.theta[i]

    def __set_view_angle_mpl(self, ax):
        phi_den_max, theta_den_max = self.__max_density_angles()
        elevation = -theta_den_max * 180.0 / np.pi
        rotation = phi_den_max * 180.0 / np.pi
        ax.view_init(elev=elevation, azim=rotation)

    def __set_view_angle_mayavi(self):
        phi_den_max, theta_den_max = self.__max_density_angles()
        if theta_den_max < np.pi / 2:
            elevation = theta_den_max * 180.0 / np.pi + 90
        else:
            elevation = theta_den_max * 180.0 / np.pi - 90
        rotation = phi_den_max * 180.0 / np.pi
        mlab.view(elevation=elevation, azimuth=rotation)

    def __rgba_colors(self, main_color, abs_square, contrast_threshold=0.10):
        if main_color == "red":
            rgb = np.array([1.0, 0.15, 0.15])
            fixed_opacity = 0.7
            slope_factor = 4.5
            x0 = 0.4
        else:
            rgb = np.array([0.3, 0.3, 1.0])
            fixed_opacity = 0.2
            slope_factor = 4.0
            x0 = 0.6
        if self.__density_contrast(abs_square) < contrast_threshold:
            opacity = np.ones([self.theta.size, self.phi.size]) * fixed_opacity
        else:
            normalized = self.__normalize_density(abs_square)
            opacity = (np.tanh((normalized - x0) * slope_factor) + 1) / 2
        rgba_colors = np.zeros([self.theta.size, self.phi.size, 4])
        for i in range(self.theta.size):
            for j in range(self.phi.size):
                rgba_colors[i, j, :3] = rgb
                rgba_colors[i, j, 3] = opacity[i, j]
        return rgba_colors

    def __plot_rgba_sphere(self, ax, radius, colors):
        x, y, z = self.__cartesian_grid(radius)
        surface = ax.plot_surface(
            x,
            y,
            z,
            rstride=1,
            cstride=1,
            facecolors=colors,
            linewidth=0,
            antialiased=False,
        )
        return surface

    def __plot_cmap_sphere(self, fig, ax, den, radius=1):
        x, y, z = self.__cartesian_grid(radius)
        den_colors = self.__normalize_density(den)
        surf = ax.plot_surface(
            x,
            y,
            z,
            rstride=1,
            cstride=1,
            facecolors=mpl.cm.jet(den_colors),
            linewidth=0,
            antialiased=False,
        )
        ax.set_axis_off()
        self.__set_view_angle_mpl(ax)
        cmap_norm = mpl.colors.Normalize(vmin=den.min(), vmax=den.max())
        fig.colorbar(
            mpl.cm.ScalarMappable(cmap=mpl.cm.jet, norm=cmap_norm),
            ax=ax,
            shrink=0.4,
        )


def run(data_path, prefix, view_tool, view_mode):
    data_plotting = DensityPlot(data_path, prefix)
    if view_tool == "matplotlib":
        if view_mode == "one-sphere":
            data_plotting.show_in_one_sphere()
        else:
            data_plotting.show_in_two_spheres()
    else:
        data_plotting.show_mayavi_oriented()
    del data_plotting


if __name__ == "__main__":
    root_dir = os.path.join(
        os.path.expanduser("~"), "programs/bec_spherical_shell/output"
    )
    p = argparse.ArgumentParser()
    p.add_argument(
        "--job-prefix",
        dest="prefix",
        type=str,
        help="prefix for file names = string before the first underscore",
    )
    p.add_argument(
        "--data-path",
        dest="data_path",
        type=str,
        default=root_dir,
        help="path to output file generated after time evolution",
    )
    p.add_argument(
        "--view-tool",
        dest="view_tool",
        type=str,
        default="matplotlib",
        help="Tool to plot : `matplotlib` or `mayavi`",
    )
    p.add_argument(
        "--view-mode",
        dest="view_mode",
        type=str,
        default="one-sphere",
        help="Whether to plot color-scale densities in one or two spheres.",
    )
    args = p.parse_args()
    run(**vars(args))
