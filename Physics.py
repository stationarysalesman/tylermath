import numpy as np

def sphereHydrodynamics(V):
    """Display hydrodynamic properties of a sphere of volume V."""
    kT = 1.38e-23 * 293.15
    eta = 1.002e-3 # kg/(m*s)
    D = kT / (6 * V * eta)
    tau = 1 / (6 * D)
    return (D, tau)


def calcDiffusionConstants(rho):
    """Calculate the diffusion constants of a disk of given axial ratio, relative to that of a sphere."""
    print('rho is {}'.format(rho))
    # calculate S' from rho
    S_prime = 0.0 
    term_a = float(rho ** 2 - 1)
    term_b = float(1 - rho ** 2)
    if rho < 1:
        a  = 1 - (rho**2)
        b = 1./np.sqrt(a)
        c = np.sqrt(a) 
        S_prime = b * np.arctan(c / rho)
    elif rho > 1:
        a = (rho ** 2) - 1
        b = 1./np.sqrt(a)
        c = np.sqrt(a)
        S_prime = b * np.log(rho + c)
    else:
        print('no111!!!')
        return [1,1] 
    # Calculate the anisotropic components of rotational diffusion, i.e. D_parallel and D_perp
    f_numerator = 3 * rho * (rho - S_prime)
    f_denom = 2 * (rho ** 2 - 1) 
    D_parallel = f_numerator / f_denom
    paren = 2 * (rho ** 2) - 1
    f_numerator = 3 * rho * (paren * S_prime - rho)
    f_denom = 2 * (np.power(rho, 4) - 1)
    D_perp = f_numerator / f_denom
    return (D_parallel, D_perp)

def calcTaus(D_parallel, D_perp):
    # Calculate relatation times from anisotropic rotational diffusion coefficients
    tau_1 = 1 / (6 * D_perp)
    tau_2 = 1 / (5 * D_perp + D_parallel)
    tau_3 = 1 / (2 * D_perp + 4 * D_parallel)
    return [tau_1, tau_2, tau_3]
   

def EllipsoidRotationalProperties():
    """Calculate diffusion coefficients and rotational relaxation times for 
    ellipsoids of revolution with some range of axial ratios."""
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    data = dict() 
    total = 60 
    split = 20 
    xs = np.array([x for x in range(total)], dtype=float)
    D_par = np.zeros(len(xs))
    D_perp = np.zeros(len(xs))
    tau1s = np.zeros(len(xs))
    tau2s = np.zeros(len(xs))
    tau3s = np.zeros(len(xs))
    xs /= float(split) 
    print('{}'.format(xs))
    for i,x in enumerate(xs):
        D_par[i], D_perp[i] = calcDiffusionConstants(x)
        tau1s[i],tau2s[i],tau3s[i] = calcTaus(D_par[i], D_perp[i])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, D_par, 'b', xs, D_perp, 'r')
    ax.set_ylabel('Rotational Diffusion / D$_{sphere}$')
    ax.set_xlabel('Axial Ratio $\\rho$')
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111) 
    ax2.plot(xs, tau1s, xs, tau2s, xs, tau3s)
    ax2.set_ylabel('Rotational correlation time / $\\tau_{sphere}$')
    ax2.set_xlabel('Axial Ratio $\\rho$')
    plt.show()
