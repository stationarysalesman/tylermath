from Physics import *
import numpy as np

# D1 discs
d = 9.7e-9
r=d/2
h = 5e-9
V = np.pi * np.power(r, 2) * h 
D,tau = sphereHydrodynamics(V)
print('D1 disk parameters are D: {} and tau:{}\n'.format(D,tau))

# E3 discs
d = 11.7e-9
r=d/2
V = np.pi * np.power(r, 2) * h 
D,tau = sphereHydrodynamics(V)
print('E1 disk parameters are D: {} and tau:{}\n'.format(D,tau))


