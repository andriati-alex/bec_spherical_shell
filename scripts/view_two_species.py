import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
from numpy import pi, sqrt

mpl.rcParams["text.usetex"] = True
mpl.rcParams["figure.subplot.hspace"] = 0.0
mpl.rcParams["figure.subplot.wspace"] = 0.0


class DensityPlot:
    def __init__(self, data_path, prefix, job_id=1):
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
            file_prefix + "_speciesA_job{}_imagtime.dat".format(job_id),
            dtype=np.complex128,
        ).reshape(self.theta.size, self.phi.size)
        self.state_b = np.loadtxt(
            file_prefix + "_speciesB_job{}_imagtime.dat".format(job_id),
            dtype=np.complex128,
        ).reshape(self.theta.size, self.phi.size)
        self.abs_square_a = abs(self.state_a) ** 2
        self.abs_square_b = abs(self.state_b) ** 2

    def show_matplotlib(self):
        # create figure
        fig = plt.figure(figsize=(8, 6))
        axa = fig.add_subplot(121, projection="3d")
        axb = fig.add_subplot(122, projection="3d")
        self.__plot_cmap_sphere(fig, axa, self.abs_square_a)
        self.__plot_cmap_sphere(fig, axb, self.abs_square_b)
        axa.set_title("$|\Psi_1(\\theta,\phi)|^2$")
        axb.set_title("$|\Psi_2(\\theta,\phi)|^2$")
        axa.set_axis_off()
        axb.set_axis_off()
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

    def __mayavi_sphere(self, den, radius=1, offset=(0, 0, 0), cmap="hot"):
        x, y, z = self.__cartesian_grid(radius, offset)
        if self.__density_contrast(den) < 0.1:
            mlab.mesh(x, y, z, color=(1, 0, 0))
        else:
            mlab.mesh(x, y, z, scalars=den, colormap=cmap)


def run(data_path, prefix, view_tool, oriented):
    data_plotting = DensityPlot(data_path, prefix)
    if view_tool not in ["matplotlib", "mayavi"]:
        raise IOError("Unrecognize visualization tool {}".format(view_tool))
    if view_tool == "matplotlib":
        data_plotting.show_matplotlib()
    else:
        if oriented:
            data_plotting.show_mayavi_oriented()
        else:
            data_plotting.show_mayavi_raw()
    del data_plotting


if __name__ == "__main__":
    default_output_dir = os.path.join(
        os.path.expanduser("~"), "programs/bec_spherical_shell/output"
    )
    p = argparse.ArgumentParser(
        usage="python %(prog)s file_name_prefix [optional_args] ",
        description="Visualize density plots in a sphere",
    )
    p.add_argument(
        "prefix",
        metavar="file_name_prefix",
        type=str,
        help="prefix for file names = string before the first underscore",
    )
    p.add_argument(
        "--data-path",
        dest="data_path",
        type=str,
        default=default_output_dir,
        help="path to output file generated after time evolution",
    )
    p.add_argument(
        "--view-tool",
        dest="view_tool",
        type=str,
        default="mayavi",
        help="Tool to plot : `matplotlib` or `mayavi` (default)",
    )
    p.add_argument(
        "-oriented",
        action="store_true",
        help="Whether to rotate view angle to align spheres",
    )
    args = p.parse_args()
    run(**vars(args))
