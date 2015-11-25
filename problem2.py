import sys
import numpy as np
import matplotlib.pyplot as plt

def gaussian(r):
    return (1/np.sqrt(2*np.pi))*np.exp(-r/2)

def solver_FE(I, radius, Nr, constant, timescale):
    # Mesh points in space
    r, dr = np.linspace(0, radius, Nr+1, retstep=True)
    # Mesh points in time
    dt = 0.5
    Nt = int(round(timescale/float(dt)))
    t = np.linspace(0, timescale, Nt+1)

    surfd = np.zeros(Nr+1)
    surfd_1 = np.zeros(Nr+1)

    # Set up initial condition
    for i in range(0, Nr+1):
        surfd_1[i] = I(r[i])

    for j in range(0, Nt):
        for n in range(1, Nr):
            surfd[n] = surfd_1[n] + (constant*dt/r[n])*(((r[n+1]**0.5 - \
            r[n]**0.5)/dr) *(surfd_1[n+1]*r[n+1]**2 - 2*(surfd_1[n]*r[n]**1) \
            + surfd_1[n-1]*r[n-1]**2))

    # Insert boundary conditions
    surfd[0] = 0.087
    surfd[Nr] = 0
    # Switch variables before next step
    surfd_1, surfd = surfd, surfd_1

    return surfd_1, r, t

if __name__ == '__main__':

    radius = 3.4
    Nr = 100
    constant = 0.3
    timescale = 1e3

    surfd, r, t = solver_FE(gaussian, radius, Nr, constant, timescale)
    print len(surfd), len(r), len(t)

    plt.plot(r, surfd, lw=2)
    plt.xlabel('Radius')
    plt.ylabel('Surface Density')
    plt.savefig('surface_density_evolution.pdf')
