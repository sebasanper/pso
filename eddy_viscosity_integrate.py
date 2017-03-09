__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'

from math import sqrt, exp, cos, pi
import time
from integrate_simpson2D import simpson_integrate2D

def ainslie(ct, u0, distance_parallel, distance_perpendicular):

    # centreline = open('centreline.dat', 'w')
    # velocity = open('velocity.dat', 'w')

    h = 0.2
    L = distance_parallel
    n = int(L / h) + 1
    Uc1 = [0.0 for x in range(n)]
    d1 = [0.0 for x in range(n)]
    I0 = 8.0  # Ambient turbulence intensity according to vanlucanee. 8% on average
    k = 0.41  # von Karman constant
    Ct = ct  # Thrust coefficient
    U0 = u0
    # dr = 0.1
    Y = distance_perpendicular
    # m = int(Y / dr)

    Dmi = Ct - 0.05 - (16.0 * Ct - 0.5) * I0 / 1000.0

    def b(deficit):  # Wake width measure
        return (3.56 * Ct / (8.0 * deficit * (1.0 - 0.5 * deficit))) ** 0.5

    def F(x):  # Factor for near and far wake
        if x >= 5.5:
            return 1.0
        if x < 5.5:
            if x >= 4.5:
                return 0.65 + ((x - 4.5) / 23.32) ** (1.0 / 3.0)
            else:
                return 0.65 - ((-x + 4.5) / 23.32) ** (1.0 / 3.0)

    def E(x1, Uf, Ud, Dm):  # Eddy viscosity term
        return F(x1) * (0.015 * b(Dm) * (Uf - Ud)) + (k ** 2.0) * I0 / 100.0

    Uc1[0] = U0 * (1.0 - Dmi)  # Boundary condition at x = 2.0
    d1[0] = Dmi
    for i in range(1, n):  # For all positions in the wake centreline direction. Recursive. Whole grid
        Uc1[i] = Uc1[i - 1] + (h * 16.0 * E(i * h, U0, Uc1[i - 1], d1[i - 1]) * (Uc1[i - 1] ** 3.0 - U0 * Uc1[i - 1] ** 2.0 - Uc1[i - 1] * U0 ** 2.0 + U0 ** 3.0) / (Uc1[i - 1] * Ct * U0 ** 2.0))
        d1[i] = 1.0 - Uc1[i] / U0

    ########### Code to calculate wake deficit at a specific point instead of the whole grid. Namely, the rotor's centrepoint.

    # U = U0 * (1.0 - d1[n-1] * exp(- 3.56 * (Y / b(d1[n-1])) ** 2.0))

    ##### Code to calculate average wake deficit in all area of the rotor ###############

    ## Define function to integrate.


    def G(r, theta):
        z = sqrt(Y ** 2.0 + r ** 2.0 + 2.0 * Y * r * cos(theta))
        gauss = U0 * (1.0 - d1[n - 1] * exp(- 3.56 * (z / b(d1[n - 1])) ** 2.0))
        return r * (U0 - gauss) ** 2.0

    A = pi * 0.5 ** 2.0  ## Unitary diameter in this program.
    U = U0 - sqrt((1.0 / A) * simpson_integrate2D(G, 0.0, 0.5, 5, 0.0, 2.0 * pi, 11))

    return 1.0 - U / U0

    # centreline.close()
    # velocity.close()

if __name__ == '__main__':
    start_time = time.time()
    print ainslie(0.811089300906, 6.38872201856, 6.99893386609, 0.122166845061)
    print("--- %s seconds ---" % (time.time() - start_time))
