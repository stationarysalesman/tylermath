import numpy as np

def calcRelaxationTimes(diam, height):
    """Calculate the relaxation times of a disk of given diameter and height."""
    
    # Some definitions    
    longitudinalAxis = height / 2.
    radialAxis = diam / 2.
    rho = longitudinalAxis / radialAxis
    T = 293.15 # temperature in Kelvin

    # Calculate D, the rotational diffusion constant of a sphere with volume of disk
    # (we will approximate volume of disk using cylinder)
    V = np.pi * np.power(radialAxis, 2) * height
    eta = 1.002e-3 # dynamic, or shear, viscosity of water, in kg/(m*s)
    kT = 1.38e-23 * 293.15
    D = kT / (6 * V * eta)
    tau_D = 1 / (6 * D)
    

    # calculate S' from rho
    S_prime = None
    term_a = float(np.power(rho, 2) - 1)
    term_b = float(1 - np.power(rho, 2))
    if rho < 1:
        term_c = np.power(term_b, -0.5)
        S_prime = term_c * np.arctan(term_c / rho)
    elif rho > 1:
        term_c = np.power(term_a, -0.5)
        term_d = np.power(term_a, 0.5)
        S_prime = term_c * np.log(rho + term_d)
    else:
        return [tau_D, tau_D, tau_D] 
    
    # Calculate the anisotropic components of rotational diffusion, i.e. D_parallel and D_perp
    f_numerator = rho * (rho - S_prime)
    f_denom = term_a
    D_parallel = (3/2.) * (f_numerator / f_denom) * D
    paren = 2 * np.power(rho, 2) - 1
    f_numerator = rho * (paren * S_prime - rho)
    f_denom = np.power(rho, 4) - 1
    D_perp = (3/2.) * (f_numerator / f_denom) * D

    # Calculate relatation times from anisotropic rotational diffusion coefficients
    tau_1 = 1 / (6 * D_perp)
    tau_2 = 1 / (5 * D_perp + D_parallel)
    tau_3 = 1 / (2 * D_perp + 4 * D_parallel)

    return [tau_1, tau_2, tau_3]
   

def testPlot():
    # test function pls ignore
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    data = dict() 
    for x in range(1, 11):
        for y in range(1, 11):
            x_nm = x * 1e-9
            y_nm = y * 1e-9
            data[(x_nm,y_nm)] = calcRelaxationTimes(x_nm, y_nm) 
    xs = np.zeros(100)
    ys = np.zeros(100)
    zs = np.zeros((100, 3))
    for i, element in enumerate(data.items()):
        t, z = element
        xs[i] = t[0]
        ys[i] = t[1]
        zs[i] = z
    tau_1s, tau_2s, tau_3s = np.hsplit(zs, 3)
    tau_1s = tau_1s.ravel() 
    tau_2s = tau_2s.ravel() 
    tau_3s = tau_3s.ravel() 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, tau_1s)
    ax.set_xlabel('Diameter')
    ax.set_ylabel('Height')
    ax.set_zlabel('Tau')
    plt.show()
