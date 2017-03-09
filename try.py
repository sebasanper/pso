__author__ = 'sebasanper'

from numpy.random import random
from math import sqrt
from wake import distance
from numpy import array
from joblib import Parallel, delayed
from jensen_simple import Jensen as fit
import time
import sys


## Inertia weight by Inertia weight strategies in particle swarm optimization by Bansal et al.
def pso_horns():
    data = open('3horns_pso.dat', 'w')
    data2 = open('3swarm_horns.dat', 'w')
    data3 = open('3globalfit.dat', 'w')

    np = 20
    nt = 45
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
        xt = 400.0 * l
        yt = 3600 * k
        return xt, yt

    ## Produce starting positions. Includes boundaries of Horns Rev
    for n in range(np):
        for tur in range(nt):
            particles[n][tur] = create()

    best_layout = particles[0]

    for i in range(nt):
        data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
    data.write('\n')

    for iter in range(2000):

        start_time2 = time.time()

        fitness = Parallel(n_jobs=8)(delayed(fit)(particles[i], nt, diam / 2.0) for i in range(np))
        # for p in range(np):
        #     fitness[p] = fit(particles[p], nt, diam / 2.0)
        for p in range(np):
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

        for p in range(np):
            ## Solving Constrained Nonlinear Optimization Problems with Particle Swarm Optimization by Xiaohui Hu and Russell Eberhart. For 1.49445 learning coefficients.
            vel[p] = (0.5 + random() / 2.0) * vel[p] + 2.0 * random() * (best_local[p] - particles[p]) + 2.0 * random() * (best_layout - particles[p])
            for t in range(nt):
                if vel[p][t][0] > 400.0:
                    vel[p][t][0] = 400.0
                if vel[p][t][0] < - 400.0:
                    vel[p][t][0] = - 400.0
                if vel[p][t][1] > 3600.0:
                    vel[p][t][1] = 3600.0
                if vel[p][t][1] < - 3600.0:
                    vel[p][t][1] = - 3600.0
            particles[p] = particles[p] + vel[p]
            for tur in range(nt):
                if particles[p][tur][1] > 3600.0:
                    particles[p][tur][1] = random() * 3600.0
                elif particles[p][tur][1] < 0.0:
                    particles[p][tur][1] = random() * 3600.0
                if particles[p][tur][0] > 400.0:
                    particles[p][tur][0] = random() * 400.0
                elif particles[p][tur][0] < 0.0:
                    particles[p][tur][0] = random() * 400.0

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

        print best_global_fitness
        if iter % 50 == 0:
            sys.stdout.flush()

        print("Iteration --- %s seconds ---" % (time.time() - start_time2))
    data.close()
    data2.close()
    data3.close()

if __name__ == '__main__':
    start_time = time.time()
    pso_horns()
    print("Total job = --- %s minutes ---" % ((time.time() - start_time) / 60.0))
