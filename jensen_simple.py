__author__ = 'sebasanper'

def Jensen(a, Nt, rad):
    from math import sqrt, log
    import wake
    # from power_curve5MW import power5MW_kW as power  # Power curve 5 MW NREL

    windrose = open('horns_rev_windrose2.dat', 'r')
    nt = Nt
    layout_x = [0.0 for x in range(nt)]
    layout_y = [0.0 for x in range(nt)]
    for x in range(nt):
        layout_x[x] = float(a[x][0])
        layout_y[x] = float(a[x][1])

    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    windrose.close()
    summation = 0.0

    # Power curve 2 MW Siemens.
    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return 0.0003234808 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 **5.0 - 30.3162345004 * U0 **4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
        else:
            return 0.0

    # for wind in range(0, len(windrose_speed)):
    for wind in range(0, 1):
        # U1 = windrose_speed[wind]  # Free stream wind speed
        # U0 = U1 * (90.0 / 10.0) ** 0.11  # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        U0 = 8.5
        k = 0.04  # Decay constant
        r0 = rad  # Turbine rotor radius
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        for turbine in range(nt):
            flag = [False for x in range(nt)]
            proportion = [0.0 for x in range(nt)]
            for i in range(nt):
                try:
                    proportion[i], flag[i] = wake.determine_if_in_wake(layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i], k, r0, angle3)
                except TypeError:
                    print layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i], k, r0, angle3

            # Matrix with effect of each turbine <i = turbine> on every other turbine <j> of the farm
            for j in range(nt):
                if turbine != j and flag[j] is True:
                    wake_deficit_matrix[turbine][j] = proportion[j] * wake.wake_deficit(0.81, k, wake.distance(layout_x[turbine], layout_y[turbine], layout_x[j], layout_y[j]), r0)
                elif turbine == j:
                    wake_deficit_matrix[turbine][j] = 0.0

        total_deficit = [0.0 for x in range(nt)]
        total_speed = [0.0 for x in range(nt)]
        for j in range(nt):
            for i in range(nt):
                total_deficit[j] += wake_deficit_matrix[i][j] ** 2.0
            total_deficit[j] = sqrt(total_deficit[j])
            total_speed[j] = U0 * (1.0 - total_deficit[j])

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        efficiency = 0.0
        for l in range(nt):
            profit += power(total_speed[l])
        efficiency = profit * 100.0 / (float(nt) * power(U0))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        summation += efficiency_proportion[wind]
    return efficiency  # Farm efficiency

# if __name__ == '__main__':
#     start_time = time.time()
#     Jensen()
# print("--- %s seconds ---" % (time.time() - start_time))