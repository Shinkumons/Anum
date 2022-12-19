from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import math

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

def plotODE(t, x1, x2):
    plt.plot(t, x1, label="m1")
    plt.plot(t, x2, label="m2")
    plt.legend(loc='best')
    plt.grid()
    plt.show()

def relative_error(ts_dx, te_dx, a, b):
    return abs(ts_dx - te_dx) <= 1e-12 * (abs(a) + abs(b)) + 1e-13

def getMaxRec(t_start, t_end, init, m1, m2, k1, k2, omega,
              ts_dx=None, te_dx=None, ts_val = None, a=0, b=100):
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




def plotDiagrams(initial_states, m1, m2, k1, k2, omega):
    t = np.linspace(0, 100, 1000)
    f = lambda x: np.sin(x*omega)
    for i, init in enumerate(initial_states):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        f_args = (m1, m2, k1, k2, f)
        sol = odeint(func, init, t, args = f_args)
        ax1.plot(t, sol[:, 0])
        ax1.plot(t, sol[:, 1])
        ax1.grid()
        ax2.plot(sol[:, 0], sol[:, 1])
        plt.grid()
        plt.show()

def plot_array_ODE(initial_states, m1, m2, k1, k2, omega):
    for init in initial_states:
        t, x1, x2 = computeODE(init, m1, m2, k1, k2, omega)
        plotODE(t, x1, x2)

def q6_compute(init, m1, m2, k1, k2, omega):
    t, x1, x2 = computeODE(init, m1, m2, k1, k2, omega)
    return max(x1)
                        
def q6(init, m1, m2, k1, k2):
    om = np.linspace(0, 2.5, 1000)
    x = []
    for i in om:
        x.append(getMaxRec(0, 100, (0, 0, 0, 0), 1, 1, 1, 1, i))
    plt.plot(om, x)
    plt.grid()
    plt.show()


initial_states = [
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (0, 0, 0, 1.1)
]
omega = 1
# plotDiagrams(initial_states, 1, 1, 1, 1, omega)

computeODE((0, 0, 0, 0), 1, 1, 1, 1, 0.62)

q6((0, 0, 0, 0), 1, 1, 1, 1)
q6((0, 0, 0, 0), 3, 1, 1, 2)


