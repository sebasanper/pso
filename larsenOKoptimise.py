__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Larsen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
# Vector a must be a wind farm layout. So size 80 x 2.
def larsen(a, windrose_angle, windrose_speed, windrose_frequency):

    import wake_geometry as wake
    from math import sqrt, log, tan, cos, pi
    from numpy import deg2rad
    nt = len(a)  # Number of turbines

    summation = 0.0

    layout_x = [0.0 for x in range(nt)]
    layout_y = [0.0 for x in range(nt)]
    for x in range(nt):
        layout_x[x] = float(a[x][0])
        layout_y[x] = float(a[x][1])

    def Ct(U0):
        if U0 < 4.0:
            return 0.1
        elif U0 <= 25.0:
            return 0.00000073139922126945 * U0 ** 6.0 - 0.0000668905596915255 * U0 ** 5.0 + 0.0023937885 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
        else:
            return 0.0

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return 0.0003234808 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 **5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
        else:
            return 0.0

    # for U0 in range(4, 20):
    r0 = 40.0  # Turbine rotor radius
    D = 2.0 * r0
    A = pi * r0 ** 2.0
    H = 70.0  # Hub height
    ia = 0.08  # Ambient turbulence intensity according to vanlucanee. 8% on average

    def deff(U0):
        return D * sqrt((1.0 + sqrt(1.0 - Ct(U0))) / (2.0 * sqrt(1.0 - Ct(U0))))

    rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
    r95 = 0.5 * (rnb + min(H, rnb))

    def x0(U0):
        return 9.5 * D / ((2.0 * r95 / deff(U0)) ** 3.0 - 1.0)

    def c1(U0):
        return (deff(U0) / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (Ct(U0) * A * x0(U0)) ** (- 5.0 / 6.0)  # Prandtl mixing length

    def distance_to_front(x, y, theta, r):
        theta = deg2rad(theta)
        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

    for wind in range(0, len(windrose_angle)):
        U1 = windrose_speed[wind]  # Free stream wind speed
        U0 = U1 * (70.0 / 10.0) ** 0.11 # Power or log law for wind shear profile
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        distance = [[0.0 for x in range(2)] for x in range(nt)]

        U = [U0 for x in range(nt)]
        total_deficit = [0.0 for x in range(nt)]

        for tur in range(nt):
            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle, 100000000.0), tur]
        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
            flag = [False for x in range(nt)]
            proportion = [0.0 for x in range(nt)]
            perpendicular_distance = [0.0 for x in range(nt)]
            parallel_distance = [0.0 for x in range(nt)]
            for i in range(turbine + 1, nt):
                proportion[distance[i][1]], flag[distance[i][1]], perpendicular_distance[distance[i][1]], parallel_distance[distance[i][1]] = wake.determine_if_in_wake(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]], A, c1(U[distance[turbine][1]]), Ct(U[distance[turbine][1]]), angle3, r0)
                if parallel_distance[distance[i][1]] > 0.0:
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[i][1]] * wake.wake_deficit(U[distance[turbine][1]], Ct(U[distance[turbine][1]]), A, parallel_distance[distance[i][1]], perpendicular_distance[distance[i][1]], c1(U[distance[turbine][1]]))
                else:
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        for l in range(nt):
            profit += power(U[l])
        efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        summation += efficiency_proportion[wind]
    return 100.0 - summation
