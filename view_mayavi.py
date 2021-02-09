"""
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
"""

import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pathlib import Path
import numpy as np
from mayavi import mlab

out_path = str(Path.home()) + "/programs/bec_spherical_shell/output/"
psi_a = np.loadtxt(
    out_path + "test2_final_speciesA_imagtime.dat", dtype=np.complex128
).reshape(201, 201)
psi_b = np.loadtxt(
    out_path + "test2_final_speciesB_imagtime.dat", dtype=np.complex128
).reshape(201, 201)

# set data
phi_values = np.linspace(0, 2 * np.pi, 201)
tht_values = np.linspace(0, np.pi, 201)
phi_grid, tht_grid = np.meshgrid(phi_values, tht_values)
r = abs(psi_a) ** 2
x = np.sin(tht_grid) * np.cos(phi_grid)
y = np.sin(tht_grid) * np.sin(phi_grid)
z = np.cos(tht_grid)

#fcolors = (r - r.min()) / (r.max() - r.min())
fcolors = (np.tanh(((r - r.min()) / (r.max() - r.min()) - 0.4) * 4) + 1) / 2
rgba_colors = np.zeros([201, 201, 4])
for i in range(201):
    for j in range(201):
        rgba_colors[i, j, 0] = 1.0
        rgba_colors[i, j, 3] = fcolors[i, j]

mlab.figure(bgcolor=(1,1,1))
mlab.mesh(x,y,z,scalars=fcolors,colormap="magma")
x = x + 3.0
mlab.mesh(x,y,z,scalars=fcolors,colormap="magma")
mlab.show()

## ax = fig.add_subplot(121, projection="3d")
#ax = fig.gca(projection="3d")
#
## Set the aspect ratio to 1 so our sphere looks spherical
## fig = plt.figure(figsize=plt.figaspect(1.))
#
#surf_a = ax.plot_surface(x, y, z, facecolors=rgba_colors, linewidth=0)
#ax.set_axis_off()
#stride = r.argmax()
#row, col = int(stride / r.shape[0]), stride % r.shape[1]
#print(tht_values[row] * 180.0 / np.pi)
#if tht_values[row] < np.pi / 2:
#    elevation = -tht_values[row] * 180.0 / np.pi
#    rotation = phi_values[col] * 180.0 / np.pi
#else:
#    elevation = tht_values[row] * 180.0 / np.pi - 90.0
#    rotation = phi_values[col] * 180.0 / np.pi
#ax.view_init(elev=elevation, azim=rotation)
#
## set data
#phi_values = np.linspace(0, 2 * np.pi, 201)
#tht_values = np.linspace(0, np.pi, 201)
#phi_grid, tht_grid = np.meshgrid(phi_values, tht_values)
#r = abs(psi_b) ** 2
#x = np.sin(tht_grid) * np.cos(phi_grid)
#y = np.sin(tht_grid) * np.sin(phi_grid)
#z = np.cos(tht_grid)
#
#fcolors = (r - r.min()) / (r.max() - r.min())
#rgba_colors = np.zeros([201, 201, 4])
#for i in range(201):
#    for j in range(201):
#        rgba_colors[i, j, 2] = 1.0
#        rgba_colors[i, j, 3] = fcolors[i, j]
#
## fig = plt.figure()
## ax = fig.add_subplot(122, projection="3d")
#
## Set the aspect ratio to 1 so our sphere looks spherical
## fig = plt.figure(figsize=plt.figaspect(1.))
#
## surf_b = ax.plot_surface(x, y, z, facecolors=cm.hot(fcolors))
#surf_b = ax.plot_surface(x, y, z, facecolors=rgba_colors)
## ax.set_axis_off()
## ax.view_init(elev=elevation, azim=rotation)
#
## Plot the surface.
## surf = ax.plot_surface(X, Y, Z, cmap="jet", linewidth=0, antialiased=False)
#
## Customize the z axis.
## ax.set_zlim(-1.01, 1.01)
## ax.zaxis.set_major_locator(LinearLocator(10))
## ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
## Add a color bar which maps values to colors.
## fig.colorbar(surf)
#
#plt.show()
