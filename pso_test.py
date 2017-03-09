__author__ = 'sebasanper'

from numpy.random import normal, random
from math import sqrt
from numpy import array

def square():
    iterations =open('iterations.dat', 'w', 1)

    np = 25
    nt = 80
    particles = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])
    vel = array([[[random(), random()] for x in range(nt)] for x in range(np)])
    fitness = [0.0 for x in range(np)]
    best_local = array([[[0.0, 0.0] for x in range(nt)] for x in range(np)])
    best_own_fitness = [0.0 for x in range(np)]
    best_global_fitness = 0.0
    for n in range(np):
        particles[n] = array([[random(), random()] for x in range(nt)])

    for iter in range(10000):
        for p in range(np):
            for t in range(nt):
                for u in range(nt):
                    fitness[p] += (particles[p][t][0] - particles[p][u][0]) ** 2.0 + (particles[p][t][1] - particles[p][u][1]) ** 2.0
                fitness[p] = sqrt(fitness[p])
            if fitness[p] > best_own_fitness[p]:
                best_own_fitness[p] = fitness[p] * 1.0
                best_local[p] = particles[p] * 1.0
            if fitness[p] > best_global_fitness:
                best_global_fitness = fitness[p] * 1.0
                best = p * 1
                for t in range(nt):
                    for j in range(2):
                        iterations.write('{0:f} '.format(particles[best][t][j]))
                    iterations.write('\n')
                iterations.write('\n')
        for p in range(np):
            for t in range(nt):
                for coord in range(2):
                    vel[p][t][coord] = 0.72984 * vel[p][t][coord] + 1.49617 * random() * (best_local[p][t][coord] - particles[p][t][coord]) + 1.49617 * random() * (particles[best][t][coord] - particles[p][t][coord])
                    particles[p][t][coord] = particles[p][t][coord] + vel[p][t][coord]
                    while particles[p][t][coord] > 1.0 or particles[p][t][coord] < 0.0:
                            particles[p][t][coord] = random()

        if iter%100==0:
            print best_global_fitness
    iterations.close()

    print particles[best]

def rosenbrock():
    rosenb = open('rosenbrock.dat', 'w')

    np = 20
    particles = array([[0.0, 0.0] for x in range(np)])
    vel = array([[0.0, 0.0] for x in range(np)])
    fitness = [0.0 for x in range(np)]
    best_local = array([[0.0, 0.0] for x in range(np)])
    for n in range(np):
        if random() < 0.5:
            sign1 = 1.0
        else:
            sign1 = - 1.0
        if random() < 0.5:
            sign2 = 1.0
        else:
            sign2 = - 1.0
        particles[n] = array([sign1 * 5.0 * random(), sign2 * 5.0 * random()])
        if random() < 0.5:
            sign3 = 1.0
        else:
            sign3 = - 1.0
        if random() < 0.5:
            sign4 = 1.0
        else:
            sign4 = - 1.0
        vel[n] = array([[5.0 * sign3 * random(), 5.0 * sign4 * random()]])
    for iter in range(2000):
        for p in range(np):
            fitness[p] = (1.0 - particles[p][0]) ** 2.0 + 100.0 * (particles[p][1] - particles[p][0] ** 2.0) ** 2.0
            if iter == 0:
                best_own_fitness = [1000.0 for x in range(np)]
                best_global_fitness = 1000.0 # fitness[0]
                best = 1
            if fitness[p] < best_own_fitness[p]:
                best_own_fitness[p] = fitness[p]
                best_local[p] = particles[p]
            if fitness[p] < best_global_fitness:
                best_global_fitness = fitness[p]
                best = p
        for p in range(np):
            vel[p] = (0.5 + random() / 2.0) * vel[p] + 2.0 * random() * (best_local[p] - particles[p]) + 2.0 * random() * (particles[best] - particles[p])
            particles[p] = particles[p] + vel[p]
            if particles[p][0] > 5.0:
                particles[p][0] = 5.0
            if particles[p][0] < - 5.0:
                particles[p][0] = - 5.0
            if particles[p][1] > 5.0:
                particles[p][1] = 5.0
            if particles[p][1] < - 5.0:
                particles[p][1] = - 5.0
        for n in range(np):
            rosenb.write('{0:f} {1:f}\n'.format(particles[n][0], particles[n][1]))
        rosenb.write('\n')
    rosenb.close()

    print best_global_fitness
    print particles[best]

if __name__ == '__main__':
    square()