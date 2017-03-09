__author__ = 'sebasanper'

import area
from math import sqrt, sin, cos, tan
from numpy import deg2rad

# Basic two turbines one behind the other . but program generalised to any crosswind distance. Geometry only

def wake_radius(r0, k, x):
    return r0 + k * x

def wake_speed(U0, Ct, k, x, r0):
    return U0 * (1.0 - (1.0 - sqrt(1.0 - Ct))/(1.0 + (k * x)/r0) ** 2.0)

def wake_deficit(Ct, k, x, r0):
    return (1.0 - sqrt(1.0 - Ct))/(1.0 + (k * x)/r0) ** 2.0

def partial_wake_deficit(Ct, k, x, r0, d):
    shadow_area = area.AreaReal(r0, wake_radius(r0, k, x), d)
    return shadow_area.area() * (1.0 - (1.0 - sqrt(1.0 - Ct))/(1.0 + (k * x)/r0))

def root_square_sum(*args):
    b = 0.0
    for i in args:
        b += args[i-1] ** 2.0
    return sqrt(b)

def distance(t1_x, t1_y, t2_x, t2_y):
    return abs(sqrt((t1_x-t2_x) ** 2.0 + (t1_y - t2_y) ** 2.0))

def crosswind_distance(w_d, t1_x, t1_y, t2_x, t2_y):
    return abs(t1_x * (- sin(w_d)) + t1_y * cos(w_d) - t2_x * (- sin(w_d)) - t2_y * cos(w_d))

def determine_if_in_wake(xt, yt, xw, yw, k, r0, alpha):  # According to Jensen Model only
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
    radius = wake_radius(r0, k, distance_to_turbine)
    # print radius
    if (xw - xt) * cos(alpha) + (yw - yt) * sin(alpha) <= 0.0:
        if abs(radius) >= abs(distance_to_centre):
            if abs(radius) >= abs(distance_to_centre) + r0:
                fraction = 1.0
                value = True
                return fraction
            elif abs(radius) < abs(distance_to_centre) + r0:
                fraction = area.AreaReal(r0, radius, distance_to_centre).area()
                value = True
                return fraction
        elif abs(radius) < abs(distance_to_centre):
            if abs(radius) <= abs(distance_to_centre) - r0:
                fraction = 0.0
                value = False
                return fraction
            elif abs(radius) > abs(distance_to_centre) - r0:
                fraction = area.AreaReal(r0, radius, distance_to_centre).area()
                value = True
                return fraction
    else:
        return 0.0


if __name__ == '__main__':
    jensen = open('../Eddy Viscosity/jensen.dat', 'w')
    U0 = 10.0
    for i in range(1000):
        radius = wake_radius(40.0, 0.04, i * 0.02 * 80.0)
        for j in range(1000):
            if j * 0.005 * 80.0 <= radius:
                Uy = wake_speed(U0, 0.8, 0.04, i * 0.02 * 80.0, 40.0)
            elif j * 0.005 * 80.0 > radius:
                Uy = U0
            jensen.write('{0:f}\t{1:f}\t{2:f}\n'.format(i * 0.02 + 2.0, j * 0.005, Uy))
        jensen.write('\n')
    jensen.close()

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
