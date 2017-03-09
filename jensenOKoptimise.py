__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
# Vector a must be a wind farm layout. So size 80 x 2.
def jensen(a, windrose_angle, windrose_speed, windrose_frequency):
    import wake
    from math import sqrt, log, tan, cos
    from numpy import deg2rad

    nt = len(a)
    layout_x = [0.0 for x in range(nt)]
    layout_y = [0.0 for x in range(nt)]
    for x in range(nt):
        layout_x[x] = float(a[x][0])
        layout_y[x] = float(a[x][1])


        # Ct curves.
        # Polynomials 6th, 4th 3rd
        # - 3.10224672352816e-5 * U0 ** 4.0 + 0.0021367624 * U0 ** 3.0 - 0.0495873986 * U0 ** 2.0 + 0.3976324804 * U0 - 0.1608576035
        # 3.374593e-4 * U0 ** 3.0 - 0.0136412226 * U0 ** 2.0 + 0.1118003309 * U0 + 0.5782039288

        # Ct look-up table for linear interpolation
        # def thrust_table(v):
            #     if v == 4: return 0.82
            #     if v == 5: return 0.81
            #     if v == 6: return 0.8
            #     if v == 7: return 0.81
            #     if v == 8: return 0.81
            #     if v == 9: return 0.78
            #     if v == 10: return 0.74
            #     if v == 11: return 0.65
            #     if v == 12: return 0.57
            #     if v == 13: return 0.41
            #     if v == 14: return 0.31
            #     if v == 15: return 0.25
            #     if v == 16: return 0.2
            #     if v == 17: return 0.17
            #     if v == 18: return 0.14
            #     if v == 19: return 0.12
            #     if v == 20: return 0.1
            #     if v == 21: return 0.09
            #     if v == 22: return 0.08
            #     if v == 23: return 0.07
            #     if v == 24: return 0.06
            #     if v == 25: return 0.05

        # 2 step thrust curve:
        # 0.80 if U0 in 4 - 9 m/s
        # 0.250625 if U0 in 10 - 25 m/s


    def Ct(U0):
        if U0 < 4.0:
            return 0.1
        elif U0 <= 25.0:
            return 7.3139922126945e-7 * U0 ** 6.0 - 6.68905596915255e-5 * U0 ** 5.0 + 2.3937885e-3 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
        else:
            return 0.0

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return 3.234808e-4 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 ** 5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
        else:
            return 0.0

            # Power curve polynomials.
        # elif U0 <= 25.0
            # - 0.0110778061 * U0 ** 5.0 + 0.8986075613 * U0 ** 4.0 - 27.2165513154 * U0 ** 3.0 + 368.8877606215 * U0 ** 2.0 - 1994.1905079276 * U0 + 3712.3986113386 #  5th degree
            # - 0.5308414162 * U0 ** 3.0 - 15.4948143381 * U0 ** 2.0 + 13.1508234816 * U0  # 3rd degree

        #  Interpolation
            #table:
    # def power_table(v):
    #     if v == 4: return 66.3
    #     if v == 5: return 152
    #     if v == 6: return 280
    #     if v == 7: return 457
    #     if v == 8: return 690
    #     if v == 9: return 978
    #     if v == 10: return 1296
    #     if v == 11: return 1598
    #     if v == 12: return 1818
    #     if v == 13: return 1935
    #     if v == 14: return 1980
    #     if v == 15: return 1995
    #     if v == 16: return 1999
    #     if v == 17: return 2000
    #     if v == 18: return 2000
    #     if v == 19: return 2000
    #     if v == 20: return 2000
    #     if v == 21: return 2000
    #     if v == 22: return 2000
    #     if v == 23: return 2000
    #     if v == 24: return 2000
    #     if v == 25: return 2000
# interpolation function
    # def interpolate(minx, miny, maxx, maxy, valx):
        # return miny + (maxy - miny) * ((valx - minx) / (maxx - minx))

        # 2 step Power curve
        # 815.0333333 kw from 4 to 13 m/s, 2000 from 13-25 m/s


    # for U0 in range(4, 20):
    summation = 0.0

    def distance_to_front(x, y, theta, r):
        theta = deg2rad(theta)
        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

    # for wind in range(0, 1):
    for wind in range(0, len(windrose_angle)):
        U1 = windrose_speed[wind]  # Free stream wind speed
        U0 = U1 * (70.0 / 10.0) ** 0.11 # Power or log law for wind shear profile
        # U0 = 8.5
        k = 0.04  # Decay constant
        r0 = 40.0  # Turbine rotor radius
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        proportion = [[0.0 for x in range(nt)] for x in range(nt)]
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
            for i in range(turbine + 1, nt):
                proportion[distance[turbine][1]][distance[i][1]] = wake.determine_if_in_wake(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]], k, r0, angle3)
                deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[turbine][1]][distance[i][1]] * wake.wake_deficit(Ct(U[distance[turbine][1]]), k, wake.distance(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]]), r0)

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        for l in range(nt):
            profit += power(U[l])
        efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0

        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        summation += efficiency_proportion[wind]

    return 100.0 - summation
