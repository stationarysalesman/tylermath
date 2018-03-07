from Physics import *
import numpy as np
import matplotlib.pyplot as plt


def disc_volume(r, h):
    return np.pi * (r ** 2) * h


def equiv_sphere_radius(volume):
    return np.power(volume / ((4/3.) * np.pi), 1/3.)


# Empty D1 discs
D1_diam = 9.7e-9
#D1_diam = 10.2e-9
D1_height = 5e-9
D1_rho = D1_height / D1_diam
D1_volume = disc_volume(D1_diam/2, D1_height)
print('Sphere hydrodynamics for D1 discs:')
D_0, tau_0=sphereHydrodynamics(D1_volume)
print('Diffusion constant: {}'.format(D_0))
print('Tau: {}'.format(tau_0))
print('Calculations based on axial ratio of {}:'.format(D1_rho))
Dpar,Dperp=calcDiffusionConstants(D1_rho)
Dpar_true = Dpar * D_0
Dperp_true = Dperp * D_0
print('D parallel: {}'.format(Dpar_true))
print('D perp: {}'.format(Dperp_true))
taus=np.array(calcTaus(Dpar_true,Dperp_true))
print('Taus: {}'.format(taus))

# Empty E3 discs
E3_diam = 11.7e-9
#E3_diam = 12.7e-9
E3_height = 5e-9
E3_rho = E3_height / E3_diam
E3_volume = disc_volume(E3_diam/2, E3_height)
print('Sphere hydrodynamics for E3 discs:')
D_0, tau_0=sphereHydrodynamics(E3_volume)
print('Diffusion constant: {}'.format(D_0))
print('Tau: {}'.format(tau_0))
print('Calculations based on axial ratio of {}:'.format(E3_rho))
Dpar,Dperp=calcDiffusionConstants(E3_rho)
Dpar_true = Dpar * D_0
Dperp_true = Dperp * D_0
print('D parallel: {}'.format(Dpar_true))
print('D perp: {}'.format(Dperp_true))
taus=np.array(calcTaus(Dpar_true,Dperp_true))
print('Taus: {}'.format(taus))

"""
# D1 discs with some KRAS bound
a = 1.55e-9
c = 3.5e-9
#kras_volume = (4/3.) * np.pi * a * a * c
kras_volume = (3.1e-9 ** 2) * 7.1e-9
print('Sphere hydrodynamics for D1 discs with 4 Kras bound (one per leaflet):')
D,t = sphereHydrodynamics(D1_volume + (2 * kras_volume))
print('D: {}, tau: {}'.format(D,t))
"""

"""
# Plot tau as a function of total volume for sphere approximation
xs = np.array([x for x in range(100)], dtype=float)
for i in range(len(xs)):
    xs[i] *= 1e-25
ys = np.zeros(len(xs))
for i, x in enumerate(xs):
    d,t = sphereHydrodynamics(x)
    ys[i] = t 

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xs, ys, 'b')
ax.set_ybound([90e-9, 600e-9])
ax.set_ylabel('Tau')
ax.set_xlabel('Volume')
plt.show()
"""
