__author__ = 'sebasanper'

from numpy.random import random
from math import sqrt
from wake import distance
from numpy import array
from joblib import Parallel, delayed
from jensenOKoptimise import jensen as fit
# from ainslieOKoptimise import ainslie as fit ##Enforce 2D minimum distance between turbines.
# from larsenOKoptimise import larsen as fit

import time

## Inertia weight 0.5+rand/2.0, by: "Inertia weight strategies in particle swarm optimization" by Bansal et al.
def pso_horns():

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

    data = open('try_best_layout_jensen.dat', 'w')
    data2 = open('try_random_layout_jensen.dat', 'w')
    data3 = open('try_best_global_fitness_jensen.dat', 'w', 1)

    np = 20 ## Number of particles in swarm
    nt = 80
    diam = 80.0
    particles = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])

    if random() < 0.5:
        sign1 = 1.0
    else:
        sign1 = - 1.0
    if random() < 0.5:
        sign2 = 1.0
    else:
        sign2 = - 1.0

    vel = array([[[sign1 * random(), sign2 * random()] for x in range(nt)] for x in range(np)])
    best_local = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])
    best_own_fitness = [0.0 for x in range(np)]
    best_global_fitness = 0.0

    def create():
        k = random()
        l = random()
        xt = 5457.0 * l
        if xt <= 412.0:
            yt = k * 3907.0 + (1.0 - k) * (- 3907.0 / 412.0 * xt + 3907.0)
        elif xt <= 5040.0:
            yt = random() * 3907.0
        else:
            yt = k * (3907.0 / 417.0 * (- xt + 5457.0))
        return xt, yt

    #  Produce starting positions. Includes boundaries of Horns Rev
    for n in range(np):
        for tur in range(nt):
            particles[n][tur] = create()

    best_layout = particles[0]

    for i in range(nt):
        data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
    data.write('\n')

    # Velocity limiting to 10% to start with, for convergence, and then increase speed.s
    k = 0.1

    for iter in range(2000):

        start_time2 = time.time()
        # fitness = Parallel(n_jobs=8)(delayed(fit)(particles[i], windrose_angle, windrose_speed, windrose_frequency) for i in range(np))
        fitness = [0.0 for x in range(np)]

        # Fitness evaluation skipped if a turbine is out of boundaries. following: BrattonKennedy07 PSO.
        for p in range(np):  # Can be parallelised in the future.
            check = True
            for tur in range(nt):
                if particles[p][tur][1] > 3907.0:
                    check = False
                    break
                elif particles[p][tur][1] < 0.0:
                    check = False
                    break
                elif particles[p][tur][1] > 3907.0 / 417.0 * (- particles[p][tur][0] + 5457.0):
                    check = False
                    break
                elif particles[p][tur][1] < - 3907.0 / 412.0 * particles[p][tur][0] + 3907.0:
                    check = False
                    break
            if check:
                fitness[p] = fit(particles[p], windrose_angle, windrose_speed, windrose_frequency)
                if fitness[p] > best_own_fitness[p]:
                    best_own_fitness[p] = fitness[p]
                    best_local[p] = particles[p]
                if fitness[p] > best_global_fitness:
                    best_global_fitness = fitness[p]
                    best_layout = particles[p]

        for i in range(nt):
            data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
            data2.write('{2:d} {0:f} {1:f}\n'.format(particles[10][i][0], particles[10][i][1], i))
        data2.write('\n')
        data.write('\n')
        data3.write('{0:f}\n'.format(best_global_fitness))

        for p in range(np):  # Can be parallelised in the future.
            ## Solving Constrained Nonlinear Optimization Problems with Particle Swarm Optimization by Xiaohui Hu and Russell Eberhart. For 1.49445 learning coefficients.
            vel[p] = (0.5 + random() / 2.0) * vel[p] + 2.0 * random() * (best_local[p] - particles[p]) + 2.0 * random() * (best_layout - particles[p])
            for t in range(nt): # Velocity half the maximum distance between boundaries. following: Particle Swarm Optimization: A Tutorial by James Blondin PSO
                if vel[p][t][0] > k * 2728.5:
                    vel[p][t][0] = k * 2728.5
                if vel[p][t][0] < - k * 2728.5:
                    vel[p][t][0] = - k * 2728.5
                if vel[p][t][1] > k * 1953.5:
                    vel[p][t][1] = k * 1953.5
                if vel[p][t][1] < - k * 1953.5:
                    vel[p][t][1] = - k * 1953.5
            particles[p] = particles[p] + vel[p]

        # Find minimum distance between turbines, and if two are closer than 1D, then randomise one of them.
        for b in range(np):
            pp = 0
            while pp == 0:
                pp = 1
                for i in range(nt):
                    for j in range(nt):
                        if i != j and distance(particles[b][i][0], particles[b][i][1], particles[b][j][0], particles[b][j][1]) < diam:
                            particles[b][j] = create()
                            pp = 0

        print("Iteration --- %s seconds ---" % (time.time() - start_time2))
    data.close()
    data2.close()
    data3.close()

if __name__ == '__main__':
    start_time = time.time()
    pso_horns()
    print("Total job = --- %s minutes ---" % ((time.time() - start_time) / 60.0))
