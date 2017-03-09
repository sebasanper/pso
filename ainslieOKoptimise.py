__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Eddy Viscosity by Ainslie wake model applied to horns rev.
# Vector a must be a wind farm layout. So size 80 x 2.
def ainslie(a, windrose_angle, windrose_speed, windrose_frequency):

    from wake import crosswind_distance
    from eddy_viscosity_integrate import ainslie
    from math import sqrt, log, cos, sin, tan
    from numpy import deg2rad

    nt = len(a)
    D = 80.0  # Diameter
    summation = 0.0

    layout_x = [0.0 for x in range(nt)]
    layout_y = [0.0 for x in range(nt)]
    for x in range(nt):
        layout_x[x] = float(a[x][0] / D)
        layout_y[x] = float(a[x][1] / D)

    def determine_front(wind_angle, x_t1, y_t1, x_t2, y_t2):
        wind_angle = deg2rad(wind_angle)
        a = (x_t2 - x_t1) * cos(wind_angle) + (y_t2 - y_t1) * sin(wind_angle)
        if a > 0.0:
            return a
        else:
            return 0.0

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return - 0.5308414162 * U0 ** 3.0 + 15.4948143381 * U0 ** 2.0 + 13.1508234816 * U0
        else:
            return 0.0

    def Ct(U0):
        if U0 < 4.0:
            return 0.1
        elif U0 <= 25.0:
            return 7.3139922126945e-7 * U0 ** 6.0 - 6.68905596915255e-5 * U0 ** 5.0 + 2.3937885e-3 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
        else:
            return 0.0

    def distance_to_front(x, y, theta, r):
        theta = deg2rad(theta)
        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

    for wind in range(0, len(windrose_angle)):
        U1 = windrose_speed[wind]  # Free stream wind speed
        U0 = U1 * (70.0 / 10.0) ** 0.11  # Power or log law for wind shear profile
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        distance = [[0.0 for x in range(2)] for x in range(nt)]
        total_deficit = [0.0 for x in range(nt)]
        total_speed = [U0 for x in range(nt)]

        for tur in range(nt):
            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle, 1000000000.0), tur]
        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += wake_deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            total_speed[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
            parallel_distance = [0.0 for x in range(0, nt)]
            perpendicular_distance = [0.0 for x in range(0, nt)]
            for i in range(turbine + 1, nt):
                parallel_distance[distance[i][1]] = determine_front(angle3, layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]])
                perpendicular_distance[distance[i][1]] = crosswind_distance(deg2rad(angle3), layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]])
                if perpendicular_distance[distance[i][1]] <= 1.7 and parallel_distance[distance[i][1]] > 0.0: ## 1.7 gives same results as a bigger distance, many times faster.
                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = ainslie(Ct(total_speed[distance[turbine][1]]), total_speed[distance[turbine][1]], parallel_distance[distance[i][1]], perpendicular_distance[distance[i][1]])
                else:
                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(len(windrose_frequency))]
        efficiency = 0.0
        for l in range(nt):
            profit += power(total_speed[l])
            efficiency = profit * 100.0 / (nt * power(total_speed[distance[0][1]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        summation += efficiency_proportion[wind]
    return 100.0 - summation
