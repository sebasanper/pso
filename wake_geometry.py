__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'

import area
from math import sqrt, sin, cos, tan, pi
from numpy import deg2rad

# Basic two turbines one behind the other . but program generalised to any crosswind distance. Geometry only

def wake_radius(c1, ct, A, x):
    return ((35.0 / 2.0 / pi) ** (1.0 / 5.0)) * ((3.0 * c1 ** 2.0) ** (1.0 / 5.0)) * ((ct * A * x) ** (1.0 / 3.0))

def wake_speed(U0, ct, A, x, r, c1):
    return U0 * (1.0 - ((ct * A * x ** (- 2.0)) ** (1.0 / 3.0)) / 9.0 * (r ** (3.0 / 2.0) * (3.0 * c1 ** 2.0 * ct * A * x) ** (- 1.0 / 2.0) - (35.0 / 2.0 / pi) ** (3.0 / 10.0) * (3.0 * c1 ** 2.0) ** (- 1.0 / 5.0)) ** 2.0)

def wake_deficit(U0, ct, A, x, r, c1):
    return 1.0 - wake_speed(U0, ct, A, x, r, c1) / U0

def distance(t1_x, t1_y, t2_x, t2_y):
    return abs(sqrt((t1_x-t2_x) ** 2.0 + (t1_y - t2_y) ** 2.0))

def crosswind_distance(w_d, t1_x, t1_y, t2_x, t2_y):
    return abs(t1_x * (- sin(w_d)) + t1_y * cos(w_d) - t2_x * (- sin(w_d)) - t2_y * cos(w_d))

def determine_if_in_wake(xt, yt, xw, yw, A, c1, ct, alpha, r0):  # According to Larsen Model only
    # Eq. of centreline is Y = tan (d) (X - Xt) + Yt
    # Distance from point to line
    alpha = deg2rad(alpha + 180)
    distance_to_centre = abs(- tan(alpha) * xw + yw + tan(alpha) * xt - yt) / sqrt(1.0 + tan(alpha) ** 2.0)
        # print distance_to_centre
    # Coordinates of the intersection between closest path from turbine in wake to centreline.
    X_int = (xw + tan(alpha) * yw + tan(alpha) * (tan(alpha) * xt - yt)) / (tan(alpha) ** 2.0 + 1.0)
    Y_int = (- tan(alpha) * (- xw - tan(alpha) * yw) - tan(alpha) * xt + yt) / (tan(alpha) ** 2.0 + 1.0)
    # Distance from intersection point to turbine
    distance_to_turbine = sqrt((X_int - xt) ** 2.0+(Y_int - yt) ** 2.0)
    # Radius of wake at that distance
    radius = wake_radius(c1, ct, A, distance_to_turbine)
    # print radius
    if (xw - xt) * cos(alpha) + (yw - yt) * sin(alpha) <= 0.0:
        if abs(radius) >= abs(distance_to_centre):
            if abs(radius) >= abs(distance_to_centre) + r0:
                fraction = 1.0
                value = True
                return fraction, value, distance_to_centre, distance_to_turbine
            elif abs(radius) < abs(distance_to_centre) + r0:
                fraction = area.AreaReal(r0, radius, distance_to_centre).area()
                value = True
                return fraction, value, distance_to_centre, distance_to_turbine
        elif abs(radius) < abs(distance_to_centre):
            if abs(radius) <= abs(distance_to_centre) - r0:
                fraction = 0.0
                value = False
                return fraction, value, distance_to_centre, distance_to_turbine
            elif abs(radius) > abs(distance_to_centre) - r0:
                fraction = area.AreaReal(r0, radius, distance_to_centre).area()
                value = True
                return fraction, value, distance_to_centre, distance_to_turbine
    else:
        return 0.0, False, distance_to_centre, distance_to_turbine


if __name__ == '__main__':
    larsen = open('larsen.dat', 'w')
    U0 = 8.5
    r0 = 40.0  # Turbine rotor radius
    D = 2.0 * r0
    A = pi * r0 ** 2.0
    ct = 0.81
    deff = D * sqrt((1.0 + sqrt(1.0 - ct)) / (2.0 * sqrt(1.0 - ct)))
    H = 70.0  # Hub height
    ia = 0.1  # Ambient turbulence intensity
    rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
    r95 = 0.5 * (rnb + min(H, rnb))
    x0 = 9.5 * D / (((2.0 * r95 / deff) ** 3.0) - 1.0)
    # a1 = (105.0 / 2.0 / pi) ** (1.0 / 5.0) * (ct * A) ** (1.0 / 3.0)
    # b1 = (1.0 / 9.0) * 3.0 ** (- 2.0 / 5.0) * (35.0 / 2.0 / pi) ** (3.0 / 5.0) * (ct * A) ** (1.0 / 3.0)
    # def Um(x):
    #     return U0 * (1.0 - b1 * c1(x) ** (- 4.0 / 5.0) * (x0(x) + x) ** (- 2.0 / 3.0))
    # def x0(x):
    #     return ((D / 2 / a1) ** (- 3.0) * (b1 / (U0 - Um(x))) ** (3.0 / 4.0) - 1.0) ** (- 1.0) * x
    c1 = (deff / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (ct * A * x0) ** (- 5.0 / 6.0)  # Prandtl mixing length
    # def c1(x):
    #     return (D / 2.0 / a1) ** (3.0 / 2.0) * x0(x) ** (- 5.0 / 6.0)
    for i in range(1, 1500):
        radius = wake_radius(c1, ct, A, i * 1.2 + x0)
        for j in range(0, 1000):
            if j * 0.48 <= radius:
                Uy = wake_speed(U0, ct, A, i * 1.2 + x0, j * 0.48, c1)
            elif j * 0.48 > radius:
                Uy = U0
            larsen.write('{0:f}\t{1:f}\t{2:f}\n'.format(i * 1.2, j * 0.48, Uy))
        larsen.write('\n')
    larsen.close()
#     print(c1, x0, r95, rnb, deff)
#     print wake_speed(8.5, 0.81, 5026.54824574, 561.0, 0.0, c1)
#     print wake_deficit(8.5, 0.81, 5026.54824574, 561.0, 0.0, c1)
#
#
# #8.5 0.81 5026.54824574 561.0 0.0 0.145785018447 61.835158772 121.6 173.2 102.670850215

    # print 'The radius of the wake is {0:f} m'.format(radius)
    # example = partial_wake_speed(10.0, 0.91, 0.04, 17.0, 1.0, 0.0)
    # print examples
    # print wake_deficit(example, 10.0)
    #
    # casa = [1, 2, 3, 4]
    # print root_square_sum(*casa)
    #
    # print crosswind_distance(deg2rad(135), 2., 2., -3., -3.)

    # print determine_if_in_wake(424534., 6151447., 423974.0, 6151447., 0.1, 40.0, 0.0)
