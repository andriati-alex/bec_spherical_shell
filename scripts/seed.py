import sys
import numpy as np
from math import sqrt
from math import factorial
from pathlib import Path
from scipy.integrate import simps
from scipy.special import sph_harm


##  Read command line arguments
Nphi = int(sys.argv[1])
Nthe = int(sys.argv[2])
lmax = int(sys.argv[3])


# Generate grid
dphi = 2 * np.pi / (Nphi - 1)
dthe = np.pi / (Nthe - 1)
phi = np.linspace(0, 2 * np.pi, Nphi)
the = np.linspace(0, np.pi, Nthe)
phi, the = np.meshgrid(phi, the)


# the argument order in sph_harm function is m, l, phi, theta
# initiate with the first spherical harmonic (constant)  then
# add others with random weights
psi = sph_harm(0, 0, phi, the)
psi[1 : psi.shape[0] - 1, 1 : psi.shape[1] - 1] = (
   psi[1 : psi.shape[0] - 1, 1 : psi.shape[1] - 1]
   + np.random.random([psi.shape[0] - 2, psi.shape[1] - 2]) / 5
)
w = (np.random.random(int((2 * lmax + 1 + 1) * (lmax + 1) / 2)) - 0.5) / 0.2
k = 0
for l in range(1, lmax + 1):
   for m in range(-l, l + 1):
       psi = psi + w[k] * sph_harm(m, l, phi, the) / sqrt(factorial(l))
       k = k + 1


# normalize to 1. Attention to sin(theta) spherical coordinates jacobian
abs2 = abs(psi) ** 2
integral_phi = np.zeros(psi.shape[0], dtype=np.float64)
for i in range(psi.shape[0]):
    integral_phi[i] = simps(abs2[i], dx=dphi)
psi = psi / np.sqrt(simps(np.sin(the[:, 0]) * integral_phi, dx=dthe))


# Record initial state A
inp_path = str(Path.home()) + "/programs/becinashell/input/"
np.savetxt(
    inp_path + "state_speciesA_init.dat", psi.reshape(Nphi * Nthe).T, fmt="%.15E"
)

psi = sph_harm(0, 0, phi, the)
psi[1 : psi.shape[0] - 1, 1 : psi.shape[1] - 1] = (
   psi[1 : psi.shape[0] - 1, 1 : psi.shape[1] - 1]
   + np.random.random([psi.shape[0] - 2, psi.shape[1] - 2]) / 5
)
w = (np.random.random(int((2 * lmax + 1 + 1) * (lmax + 1) / 2)) - 0.5) / 0.1
k = 0
for l in range(1, lmax + 1):
   for m in range(-l, l + 1):
       psi = psi + w[k] * sph_harm(m, l, phi, the) / sqrt(factorial(l))
       k = k + 1


# normalize to 1. Attention to sin(theta) spherical coordinates jacobian
abs2 = abs(psi) ** 2
integral_phi = np.zeros(psi.shape[0], dtype=np.float64)
for i in range(psi.shape[0]):
    integral_phi[i] = simps(abs2[i], dx=dphi)
psi = psi / np.sqrt(simps(np.sin(the[:, 0]) * integral_phi, dx=dthe))


# Record initial state A
inp_path = str(Path.home()) + "/programs/becinashell/input/"
np.savetxt(
    inp_path + "state_speciesB_init.dat", psi.reshape(Nphi * Nthe).T, fmt="%.15E"
)


# Configure domain
f = open(inp_path + "state_domain.dat", "w")
f.write("{} {} {} {}".format(Nphi, Nthe, 0.002, 10000))
f.close()


# Default values are set for equation parameters.
f = open(inp_path + "state_eq.dat", "w")
f.write("-0.5 0.0 0.5 0.5 10.0 10.0 5.0")
f.close()
