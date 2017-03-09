__author__ = 'Sebastian Sanchez Perez-Moreno. Email: s.sanchezperezmoreno@tudelft.nl'

from numpy.random import random
from wake import distance
import copy
from math import sqrt
from numpy import array
from joblib import Parallel, delayed
from jensenOKoptimise import jensen as fit
#from ainslieOKoptimise import ainslie as fit
# from larsenOKoptimise import larsen as fit
import time


# def fit(particle, nt):
#     fitness = 0.0
#     for t in range(nt):
#         for u in range(nt):
#             fitness += (particle[t][0] - particle[u][0]) ** 2.0 + (particle[t][1] - particle[u][1]) ** 2.0
#         return 1.0 / sqrt(fitness)


# Inertia weight 0.5+rand/2.0, by: "Inertia weight strategies in particle swarm optimization" by Bansal et al.
def pso_horns():

    windrose = open('horns_rev_windrose.dat', 'r')
    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    windrose.close()

    data = open('ref3_layout.dat', 'w', 1)
    data2 = open('ref3_random_layout.dat', 'w', 1)
    data3 = open('ref3_best_global_fitness.dat', 'w', 1)
    allfit = open('ref3_all_fitnesses.dat', 'w', 1)
    #  pso1 first run. Np 40 (max recommended). Ainslie model combination 30 from MCDA. 20 000 iterations to see what happens. 12.91 = 87.09% eff.
    #  pso2 Less particles down to 25. Ive changed the inertia and acceleration coefficients multiplied by 0.1. added e-1. For smaller movements and less chaos. Didnt work because of e-1 i think, it moves too slow so nothing changes. 13.46 = 86.54% eff.
    #  pso3 chagned e-1 to decimal notation of the accel. coefficients. 5 particles only to test.
    #  pso4     n_iter = 80    np = 25 ## Number of particles in swarm. 2000 function calls. 8 hrs
    # pso6 same as pso 4 more iterations to  160.
    #  pso8 includes one particle as the original horns rev. 10 particles to test.
    #  ref1 reflecting, 10 particles. Changed vel and particles update with one operation per dimension. 138 seconds per iteration.
    #  ref2 10 particles, changed windrose to 30 degrees for speed. 3.5 seconds per iteration. 200 iterations.
    #  ref3 same as ref2, 2000 iterations.
    n_iter = 2000
    np = 10 ## Number of particles in swarm
    nt = 80
    diam = 80.0
    particles = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])

    # if random() < 0.5:
    #     sign1 = 1.0
    # else:
    #     sign1 = - 1.0
    # if random() < 0.5:
    #     sign2 = 1.0
    # else:
    #     sign2 = - 1.0

    vel = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])

    best_own_fitness = [100.0 for x in range(np)]
    best_global_fitness = 100.0

    def create():
        k = random()
        l = random()
        xt = 5457.0 * l
        if xt <= 412.0:
            yt = k * 3907.0 + (1.0 - k) * (- 3907.0 / 412.0 * xt + 3907.0)
        elif xt <= 5040.0:
            yt = k * 3907.0
        else:
            yt = k * (3907.0 / 417.0 * (- xt + 5457.0))
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
    # particles[np - 1] = horns
    fitness = Parallel(n_jobs=-1)(delayed(fit)(particles[i], windrose_angle, windrose_speed, windrose_frequency) for i in range(np))
    # fitness = Parallel(n_jobs=-1)(delayed(fit)(particles[i], nt) for i in range(np))
    best_layout = copy.deepcopy(particles[fitness.index(min(fitness))])
    best_local = copy.deepcopy(particles)

    for ite in range(n_iter):
        start_time2 = time.time()

        for p in range(np):  # Can be parallelised in the future.
            # Solving Constrained Nonlinear Optimization Problems with Particle Swarm Optimization by Xiaohui Hu and Russell Eberhart. For 1.49445 learning coefficients.
            for t in range(nt):
                for coord in range(2):
                    vel[p][t][coord] = 0.72984 * vel[p][t][coord] + 1.49617 * random() * (best_local[p][t][coord] - particles[p][t][coord]) + 1.49617 * random() * (best_layout[t][coord] - particles[p][t][coord])
                    particles[p][t][coord] += vel[p][t][coord]
                j = 1.0
                w = 1.0
                while particles[p][t][1] > 3907.0 or particles[p][t][1] < 0.0:
                    if particles[p][t][1] > 3907.0:
                        particles[p][t][1] = 3907.0 * 2.0 - particles[p][t][1]
                        j = random()
                    elif particles[p][t][1] < 0.0:
                        particles[p][t][1] = - particles[p][t][1]
                        j = random()
                while particles[p][t][1] < - 3907.0 / 412.0 * (particles[p][t][0]+10.) + 3907.0 or particles[p][t][1] > 3907.0 / 417.0 * (- particles[p][t][0] + 5457.0+10.):
                    if particles[p][t][1] < - 3907.0 / 412.0 * (particles[p][t][0]+10.) + 3907.0:
                        particles[p][t][0] = 2.0 * (412.0 / 3907.0) * (3907.0 - particles[p][t][1]) - particles[p][t][0]
                        w = random()
                    elif particles[p][t][1] > 3907.0 / 417.0 * (- particles[p][t][0] + 5457.0+10.):
                        particles[p][t][0] = 2.0 * (5457.0 - particles[p][t][0] * 417.0 / 3907.0) - particles[p][t][0]
                        w = random()
                vel[p][t] *= [- w, - j]
        # Find minimum distance between turbines, and if two are closer than 1D, then randomise one of them. TODO Must be 2D if used with Ainslie model, since it does not guarantee any data before 2D.
        for b in range(np):
            pp = 0
            while pp == 0:
                pp = 1
                for i in range(nt):
                    for j in range(nt):
                        if i != j and distance(particles[b][i][0], particles[b][i][1], particles[b][j][0], particles[b][j][1]) < 2.0 * diam:
                            particles[b][j] = create()
                            pp = 0


        # More 'repair' methods for particles out of boundaries shown in PhD thesis Helwig2010.
        # Chu2011 proves that the reflecting boundary method is better than random or absorbing boundary.

        for fg in range(np):
            allfit.write('{0:f}\n'.format(fitness[fg]))
        allfit.write('\n')

        fitness = Parallel(n_jobs=-1)(delayed(fit)(particles[i], windrose_angle, windrose_speed, windrose_frequency) for i in range(np))
        # fitness = Parallel(n_jobs=-1)(delayed(fit)(particles[i], nt) for i in range(np))

        for p in range(np):
            if fitness[p] < best_own_fitness[p]:
                best_own_fitness[p] = fitness[p]
                best_local[p] = copy.deepcopy(particles[p])
            if fitness[p] < best_global_fitness:
                best_global_fitness = fitness[p]
                best_layout = copy.deepcopy(particles[p])
                for i in range(nt):
                    data.write('{2:d} {0:f} {1:f}\n'.format(best_layout[i][0], best_layout[i][1], i))
                data.write('\n')
                data3.write('{1:d} {0:f}\n'.format(best_global_fitness, ite))

        for i in range(nt):
            data2.write('{2:d} {0:f} {1:f}\n'.format(particles[2][i][0], particles[2][i][1], i))
        data2.write('\n')

        print(" --- %s seconds ---" % (time.time() - start_time2))
    data.close()
    data2.close()
    data3.close()
    allfit.close()
if __name__ == '__main__':
    start_time = time.time()
    pso_horns()
    print("Full optimisation  --- %s minutes ---" % ((time.time() - start_time) / 60.0))
