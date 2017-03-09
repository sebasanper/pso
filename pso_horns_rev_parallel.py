__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
from numpy.random import random
from wake import distance
from numpy import array
from joblib import Parallel, delayed
# from jensenOKoptimise import jensen as fit
from ainslieOKoptimise import ainslie as fit
# from larsenOKoptimise import larsen as fitW
from copy import deepcopy

import time

def movement(vel, best_local, particle, best_layout, k, nt):
    for number in range(nt):
        for coord in range(2):
            vel[number][coord] = 0.72984 * vel[number][coord] + 1.49617 * random() * (best_local[number][coord] - particle[number][coord]) + 1.49617 * random() * (best_layout[number][coord] - particle[number][coord])
    for t in range(nt): # Velocity half the maximum distance between boundaries. following: Particle Swarm Optimization: A Tutorial by James Blondin PSO
        if vel[t][0] > k * 2728.5:
            vel[t][0] = k * 2728.5
        elif vel[t][0] < - k * 2728.5:
            vel[t][0] = - k * 2728.5
        if vel[t][1] > k * 1953.5:
            vel[t][1] = k * 1953.5
        elif vel[t][1] < - k * 1953.5:
            vel[t][1] = - k * 1953.5
    return particle + vel

def efficiency(particles2, windrose_angle, windrose_speed, windrose_frequency, nt):
    check = True
    for tur in range(nt):
        if particles2[tur][1] > 3907.0 or particles2[tur][1] < 0.0 or particles2[tur][1] > 3907.0 / 417.0 * (- particles2[tur][0] + 5457.0 + 10.0) or particles2[tur][1] < - 3907.0 / 412.0 * (particles2[tur][0] + 10.0) + 3907.0:
            check = False
            break
    if check:
        return fit(particles2, windrose_angle, windrose_speed, windrose_frequency)
    else:
        return 100.0

## Inertia weight 0.5+rand/2.0, by: "Inertia weight strategies in particle swarm optimization" by Bansal et al.
def pso_horns():
    time4 = time.time()
    windrose = open('horns_rev_windrose2.dat', 'r')
    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    windrose.close()

    data = open('try6_layout_jensen.dat', 'w')
    data2 = open('try6_random_layout_jensen.dat', 'w')
    data3 = open('try6_best_global_fitness_jensen.dat', 'w', 1)
    allfit = open('try6_all_fitnesses.dat', 'w', 1)

    np = 10 ## Number of particles in swarm
    nt = 80
    diam = 80.0
    particles = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])

    vel = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])
    best_own_fitness = [100.0 for x in range(np)]

    def create():
        kk = random()
        l = random()
        xt = 5457.0 * l
        if xt <= 412.0:
            yt = kk * 3907.0 + (1.0 - kk) * (- 3907.0 / 412.0 * (xt + 10.0) + 3907.0)
        elif xt <= 5040.0:
            yt = kk * 3907.0
        else:
            yt = kk * (3907.0 / 417.0 * (- xt + 5457.0 + 10.0))
        return xt, yt

    ## Produce starting positions. Includes boundaries of Horns Rev
    for n in range(np):
        for tur in range(nt):
            particles[n][tur] = create()

    # layout = open('horns_rev.dat', 'r')
    # ll = 0
    # horns = array([[0, 0] for gf in range(nt)])
    # for line in layout:
    #     columns = line.split()
    #     horns[ll] = [float(columns[0]) - 423974.0, float(columns[1]) - 6147543.0]
    #     ll += 1
    # layout.close()

    # particles[np - 1] = deepcopy(horns)

    best_local = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])
    # Fitness evaluation skipped if a turbine is out of boundaries. following: BrattonKennedy07 PSO.
    # More 'repair' methods for particles out of boundaries shown in PhD thesis Helwig2010.
    # Chu2011 proves that the reflecting boundary method is better than random or absorbing boundary. TODO implement
    fitness = Parallel(n_jobs=-1)(delayed(efficiency)(particles[i], windrose_angle, windrose_speed, windrose_frequency, nt) for i in range(np))
    print fitness
    for fg in range(np):
        allfit.write('{0:f}\n'.format(fitness[fg]))
    allfit.write('\n')

    for p in range(np):
        best_own_fitness[p] = deepcopy(fitness[p])
        best_local[p] = deepcopy(particles[p])
    best_global_fitness = min(fitness)
    print best_global_fitness
    print best_own_fitness
    best_layout = deepcopy(particles[fitness.index(min(fitness))])
    print best_local[5][54][1] - particles[5][54][1]
    # for i in range(nt):
    #     data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
    # data.write('\n')

    # Velocity limiting to 10% to start with, for convergence, and then increase speed.
    k = 1.0
    for ite in range(200):

        start_time2 = time.time()

        particles = Parallel(n_jobs=-1)(delayed(movement)(vel[i], best_local[i], particles[i], best_layout, k, nt) for i in range(np))

        # Find minimum distance between turbines, and if two are closer than 1D, then randomise one of them.
        for b in range(np):
            pp = 0
            while pp == 0:
                pp = 1
                for i in range(nt):
                    for j in range(nt):
                        if i != j and distance(particles[b][i][0], particles[b][i][1], particles[b][j][0], particles[b][j][1]) < 2.0 * diam:
                            particles[b][j] = create()
                            pp = 0

        # Fitness evaluation skipped if a turbine is out of boundaries. following: BrattonKennedy07 PSO.
        # More 'repair' methods for particles out of boundaries shown in PhD thesis Helwig2010.
        # Chu2011 proves that the reflecting boundary method is better than random or absorbing boundary. TODO implement
        fitness = Parallel(n_jobs=-1)(delayed(efficiency)(particles[i], windrose_angle, windrose_speed, windrose_frequency, nt) for i in range(np))

        for fg in range(np):
            allfit.write('{0:f}\n'.format(fitness[fg]))
        allfit.write('\n')

        for p in range(np):
            if fitness[p] < best_own_fitness[p]:
                best_own_fitness[p] = deepcopy(fitness[p])
                best_local[p] = deepcopy(particles[p])
            if fitness[p] < best_global_fitness:
                best_global_fitness = deepcopy(fitness[p])
                best_layout = deepcopy(particles[p])

        for i in range(nt):
            data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
            data2.write('{2:d} {0:f} {1:f}\n'.format(particles[1][i][0], particles[1][i][1], i))
        data2.write('\n')
        data.write('\n')
        data3.write('{0:f}\n'.format(best_global_fitness))


        print("Iteration --- %s seconds ---" % (time.time() - start_time2))
    data.close()
    data2.close()
    data3.close()
    allfit.close()

if __name__ == '__main__':
    start_time = time.time()
    pso_horns()
    print("Total job = --- %s minutes ---" % ((time.time() - start_time) / 60.0))
