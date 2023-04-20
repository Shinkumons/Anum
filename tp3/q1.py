from scipy.integrate import odeint
import numpy as np
import math
import sys

# t is always the time range
# initial_state and init are the initial value of the ode is position and
# velocity of each masses
# m1, m2, k1, k2 are the parametter of the ODE


def func(y, t, m1, m2, k1, k2, f):
    """
    Represent the ODE, should not be modified
    """
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
    Compute the ODE an return time and position of the masses
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

def relative_error(ts_dx, te_dx, a, b, eps, delt):
    # return True if the derivative between t_start and t_end are close enought
    # return False otherise
    return abs(ts_dx - te_dx) <= eps * (abs(a) + abs(b)) + delt

def getMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega):
    """
    return the max value of the local maximums returned by getLocMaxRec
    """
    locMaxi = getLocMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega)
    maxi = max([i[1] for i in locMaxi])
    return maxi

def getLocMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega,
                 ts_dx=None, te_dx=None, ts_val = None, a=0, b=100):
    """
    Compute recursively the local maximums of x1 for a given initial states and
    parametters
    t_start, t_end ->  the range of computation
    ts_dx, te_dx -> the derivative of t_start an t_end
    ts_val -> the value of t_start
    Returns a list with all locals maximums of the ODE
    """
    # base case
    if ts_val != None and relative_error(ts_dx, te_dx, a, b, 1e-12, 1e-13):
        return [(t_start, ts_val)]
    else: # recursion case
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
    """
    Compute max value of x1 for omega in [0, 2.5]
    """
    om = np.linspace(0, 2.5, 1000)
    x = []
    for i in om:
        x.append(getMaxRec(0, 100, (0, 0, 0, 0), m1, m2, k1, k2, i))
    return om, x

def getUnimodalIntervals(om, x):
    """
    Calculate approximatively strictly unimodal interval of the function
    for golden section search
    """
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

def goldenSectionSearch(a, b, m1, m2, k1, k2):
    """
    Golden section search return the maximum(resp. minimum) on a strictly
    unimodal function.
    """
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
    # Bad input management
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
