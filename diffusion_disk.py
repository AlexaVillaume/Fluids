import sys
import numpy as np
import matplotlib.pyplot as plt

def gaussian(r):
    return (1/np.sqrt(2*np.pi))*np.exp(-r/2)

def solve(initial_conds, radius, num_r, constant, timescale):
    # Grid points in space
    r, dr = np.linspace(0, radius, num_r+1, retstep=True)
    # Grid points in time
    dt = 0.1
    Nt = int(round(timescale/float(dt)))
    t = np.linspace(0, timescale, Nt+1)

    surfd = np.zeros(num_r+1)
    surfd_1 = np.zeros(num_r+1)

    # Set up initial condition
    for i in range(0, num_r+1):
        surfd_1[i] = initial_conds(r[i])

    for j in range(0, len(t)):
        for n in range(1, num_r):
            surfd[n] = surfd_1[n] + (3*constant*dt/r[n])*(((r[n+1]**0.5 - \
            r[n]**0.5)/dr) *(surfd_1[n+1]*r[n+1]**2 - 2*(surfd_1[n]*r[n]**1) \
            + surfd_1[n-1]*r[n-1]**2))

        # Insert boundary conditions
        surfd[0] = 0.087
        surfd[num_r] = 0
        # Switch variables before next step
        surfd_1, surfd = surfd, surfd_1

    return surfd_1, r, t

if __name__ == '__main__':

    alpha = 10e-2
    sound_speed = 0.6
    radius = 3.4
    num_r = 100
    constant = 0.1
    timescale = 1e3

    surfd, r, t = solve(gaussian, radius, num_r, constant, timescale)
    print len(surfd), len(r), len(t)

    plt.plot(r, surfd, lw=2)
    plt.xlabel('Radius')
    plt.ylabel('Surface Density')
    plt.savefig('surface_density_evolution.pdf')
