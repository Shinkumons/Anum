from scipy.integrate import odeint
import numpy as np
import math
import sys

def func(y, t, m1, m2, k1, k2, f):
    x1, x2, dx1, dx2 = y
    equation = [
        dx1,
        dx2,
        (1 / m1)*(-k1 * x1 + k2 * (x2 - x1)),
        (1 / m2)*(f(t) - k2 * (x2 - x1))
    ]
    return equation

def computeODE(initial_state, m1, m2, k1, k2, omega):
    """
    initial_state -> tuple of 4 element determining x1, x2, dx1/dt and dx2/dt
    k1 -> recall constant of m1
    k2 -> recall constant of m2
    m1 -> mass of the m1 block
    m2 -> mass of the m2 block
    """

    t = np.linspace(0, 100, 1000)
    f = lambda x: np.sin(x*omega)
    f_args = (m1, m2, k1, k2, f)
    sol = odeint(func, initial_state, t, args=f_args)
    return (t, sol[:, 0], sol[:, 1])

def relative_error(ts_dx, te_dx, a, b):
    return abs(ts_dx - te_dx) <= 1e-12 * (abs(a) + abs(b)) + 1e-13

def getMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega):
    locMaxi = getLocMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega)
    maxi = max([i[1] for i in locMaxi])
    return maxi

def getLocMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega,
                 ts_dx=None, te_dx=None, ts_val = None, a=0, b=100):
    if ts_val != None and relative_error(ts_dx, te_dx, a, b):
        return [(t_start, ts_val)]
    else:
        finded_max = []
        t = np.linspace(t_start, t_end, 100)
        f = lambda x: np.sin(x*omega)
        f_args = (m1, m2, k1, k2, f)
        sol = odeint(func, init, t, args=f_args)
        x1, dx1 = sol[:, 0], sol[:, 2]
        for i, val in enumerate(dx1[:-1]):
            if dx1[i] > 0 and dx1[i+1] < 0:
                init_rec = (sol[:, 0][i], sol[:, 1][i],
                            sol[:, 2][i], sol[:, 3][i])
                locmax = getLocMaxRec(t[i], t[i+1], init_rec,
                                   m1, m2, k1, k2, omega,
                                   dx1[i], dx1[i+1], x1[i])
                for j in locmax:
                    finded_max.append(j)
            elif dx1[i] == 0:
                finded_max.append((i, x1[i]))
        return finded_max

def evalOnInterval(init, m1, m2, k1, k2):
    om = np.linspace(0, 2.5, 1000)
    x = []
    for i in om:
        x.append(getMaxRec(0, 100, (0, 0, 0, 0), m1, m2, k1, k2, i))

    return om, x

def getUnimodalIntervals(om, x):
    intervals = []
    start = 0
    isInt = False
    end = 0
    for i, item in enumerate(x):
        if item >= 5 and not isInt:
            start = om[i]
            isInt = True
        if item <= 5 and isInt:
            end = om[i]
            isInt = False
            intervals.append((start, end))

    return intervals

gr = (math.sqrt(5) + 1) / 2

gr = (math.sqrt(5) + 1) / 2

def goldenSectionSearch(a, b, m1, m2, k1, k2):
    c = b - (b - a) / gr
    d = a + (b - a) / gr

    f = lambda x: getMaxRec(0, 100, (0, 0, 0, 0), m1, m2, k1, k2, x)

    while abs(b - a) > 1e-7:
        if f(c) > f(d):
            b = d
        else:
            a = c

        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (b + a) / 2


if __name__ == "__main__":
    k = len(sys.argv)
    if k > 5 :
        raise IOError(f"too many arguments, needed : 5, got : {k}")
    if k < 5 :
        raise IOError(f"Not enough argument, needed : 5, got : {k}")
    m1, m2, k1, k2 = map(float, sys.argv[1:])
    om, x = evalOnInterval((0, 0, 0, 0), m1, m2, k1, k2)
    intervals = getUnimodalIntervals(om, x)
    for start, end in intervals:
        omega = goldenSectionSearch(start, end, m1, m2, k1, k2)
        print(omega)
