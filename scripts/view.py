import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class DensityPlot:
    def __init__(self, data_path, prefix, view_tool="matplotlib"):
        file_prefix = os.path.join(data_path, prefix)
        eq_data = np.loadtxt(file_prefix + "_equation_imagtime.dat")
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
            file_prefix + "_final_speciesA_imagtime.dat", dtype=np.complex128
        ).reshape(self.theta.size, self.phi.size)
        self.state_b = np.loadtxt(
            file_prefix + "_final_speciesB_imagtime.dat", dtype=np.complex128
        ).reshape(self.theta.size, self.phi.size)
        self.abs_square_a = abs(self.state_a) ** 2
        self.abs_square_b = abs(self.state_b) ** 2
        self.view_tool = view_tool

    def show_one_sphere(self):
        phi_grid, tht_grid = np.meshgrid(self.phi, self.theta)
        # sphere of unit radius
        x = np.sin(tht_grid) * np.cos(phi_grid)
        y = np.sin(tht_grid) * np.sin(phi_grid)
        z = np.cos(tht_grid)
        # create figure
        fig = plt.figure(figsize=plt.figaspect(1.0))
        ax = fig.gca(projection="3d")
        colors_a = self.__get_colors("red", self.abs_square_a)
        colors_b = self.__get_colors("blue", self.abs_square_b)
        surf_a = ax.plot_surface(x, y, z, facecolors=colors_a, linewidth=0)
        surf_b = ax.plot_surface(x, y, z, facecolors=colors_b, linewidth=0)
        ax.set_axis_off()
        self.__set_view_angle(ax)
        plt.show()

    def __set_view_angle(self, ax):
        stride = self.abs_square_a.argmax()
        i, j = int(stride / self.theta.size), stride % self.phi.size
        if self.theta[i] < np.pi / 2:
            elevation = -self.theta[i] * 180.0 / np.pi
        else:
            elevation = self.theta[i] * 180.0 / np.pi - 90.0
        rotation = self.phi[j] * 180.0 / np.pi
        ax.view_init(elev=elevation, azim=rotation)

    def __get_colors(self, main_color, abs_square, contrast_threshold=0.10):
        max_den = abs_square.max()
        min_den = abs_square.min()
        if main_color == "red":
            color_index = 0
            fixed_opacity = 0.8
            slope_factor = 4
            x0 = 0.4
        elif main_color == "green":
            color_index = 1
            fixed_opacity = 0.6
            slope_factor = 3.5
            x0 = 0.5
        else:
            color_index = 2
            fixed_opacity = 0.4
            slope_factor = 3.0
            x0 = 0.6
        if 1.0 - min_den / max_den < contrast_threshold:
            opacity = np.ones([self.theta.size, self.phi.size]) * fixed_opacity
        else:
            normalized = (abs_square - min_den) / (max_den - min_den)
            opacity = (np.tanh((normalized - x0) * slope_factor) + 1) / 2
        rgba_colors = np.zeros([self.theta.size, self.phi.size, 4])
        for i in range(self.theta.size):
            for j in range(self.phi.size):
                rgba_colors[i, j, color_index] = 1.0
                rgba_colors[i, j, 3] = opacity[i, j]
        return rgba_colors


def run(data_path, prefix, view_tool):
    data_class = DensityPlot(data_path, prefix)
    data_class.show_one_sphere()
    del data_class


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
    args = p.parse_args()
    run(**vars(args))
