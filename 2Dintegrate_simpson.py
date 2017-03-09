__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'

#### Integrate in 2D in polar coordinates. Multiplication by r must be included in function. Function of first radius, and then angle.
def simpson_integrate2D(f, a, b, n, c, d, m):
    k = (d - c) / float(m)
    h = (b - a) / float(n)

    integral = 0.0
    for j in range(m):
        res = 0.0
        for i in range(n):
            res += abs(f(a + i * h, d + j * k) + 4.0 * (f(a + i * h + h / 2.0, d + j * k)) + f(a + (i + 1.0) * h, d + j * k))
        res *= h / 6.0
        integral += res * k
    return integral

if __name__ == '__main__':
    from math import sin, pi, sqrt
    def func(a, b):
        return 4.0 * a

    print simpson_integrate2D(func, 0.0, 1.0, 50, 0.0, 2.0 * pi, 50)
