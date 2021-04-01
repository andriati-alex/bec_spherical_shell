import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import pi


class TimeDependentState:
    def __init__(self, data_path, prefix):
        file_prefix = os.path.join(data_path, prefix)
        eq_filename = file_prefix + "_2species_equation_realtime.dat"
        if not os.path.isfile(eq_filename):
            raise IOError("\nEquation file {} not found\n".format(eq_filename))
        eq_data = np.loadtxt(file_prefix + "_2species_equation_realtime.dat")
        self.phi = np.linspace(0, 2 * pi, int(eq_data[0]))
        self.theta = np.linspace(0, pi, int(eq_data[1]))
        self.raw_state_a = np.loadtxt(
            file_prefix + "_speciesA_realtime.dat",
            dtype=np.complex128,
        )
        self.raw_state_b = np.loadtxt(
            file_prefix + "_speciesB_realtime.dat",
            dtype=np.complex128,
        )
        self.density_a = abs(self.raw_state_a) ** 2
        self.density_b = abs(self.raw_state_b) ** 2
        self.den_max = max(self.density_a.max(), self.density_b.max())
        self.den_min = min(self.density_a.min(), self.density_b.min())
        obs_data = np.loadtxt(file_prefix + "_obs_realtime.dat")
        self.time = obs_data[:, 0]
        self.overlap = obs_data[:, 2]

    def frame_min_overlap(self):
        return self.overlap.argmin()

    def frame_from_time(self, t):
        return abs(self.time - t).argmin()

    def get_density(self, frame_index=0):
        den_a = self.density_a[frame_index].reshape(self.ntheta, self.nphi)
        den_b = self.density_b[frame_index].reshape(self.ntheta, self.nphi)
        return den_a, den_b

    def preview_frame_density(self, frame_index):
        den_a = self.density_a[frame_index].reshape(
            self.theta.size, self.phi.size
        )
        den_b = self.density_b[frame_index].reshape(
            self.theta.size, self.phi.size
        )
        fig = plt.figure(figsize=(8, 6))
        axa = fig.add_subplot(121, projection="3d")
        axb = fig.add_subplot(122, projection="3d")
        self.__plot_cmap_sphere(fig, axa, den_a)
        self.__plot_cmap_sphere(fig, axb, den_b)
        self.__set_view_angle(axa, den_a)
        self.__set_view_angle(axb, den_a)
        axa.set_title("$|\Psi_1(\\theta,\phi)|^2$")
        axb.set_title("$|\Psi_2(\\theta,\phi)|^2$")
        plt.show()

    def cartesian_grid(self, radius=1, offset=(0, 0, 0)):
        phi_grid, tht_grid = np.meshgrid(self.phi, self.theta)
        # sphere of unit radius
        x = radius * np.sin(tht_grid) * np.cos(phi_grid) + offset[0]
        y = radius * np.sin(tht_grid) * np.sin(phi_grid) + offset[1]
        z = radius * np.cos(tht_grid) + offset[2]
        return (x, y, z)

    def normalize_density(self, den):
        return (den - self.den_min) / (self.den_max - self.den_min)

    def max_density_angles(self, den):
        stride = den.argmax()
        i, j = stride // self.theta.size, stride % self.phi.size
        return self.phi[j], self.theta[i]

    def __set_view_angle(self, ax, den):
        phi_den_max, theta_den_max = self.max_density_angles(den)
        elevation = -theta_den_max * 180.0 / pi
        rotation = phi_den_max * 180.0 / pi
        ax.view_init(elev=elevation, azim=rotation)

    def __plot_cmap_sphere(self, fig, ax, den, radius=1):
        x, y, z = self.cartesian_grid(radius)
        den_colors = self.normalize_density(den)
        ax.plot_surface(
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
        cmap_norm = mpl.colors.Normalize(vmin=self.den_min, vmax=self.den_max)
        fig.colorbar(
            mpl.cm.ScalarMappable(cmap=mpl.cm.jet, norm=cmap_norm),
            ax=ax,
            shrink=0.4,
        )
