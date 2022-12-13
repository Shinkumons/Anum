from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

def func(y, t, m1, m2, k1, k2, f):
    x1, x2, dx1, dx2 = y
    equation = [
        dx1,
        dx2,
        (1 / m1)*(-k1 * x1 + k2 * (x2 - x1)),
        (1 / m2)*(f(t) - k2 * (x2 - x1))
    ]
    return equation

def testODE(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt


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


initial_state_ls = [
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (1, 1, 0, 0),
    (0, 0, 1, 0),
    (1, 0, 1, 0),
    (0, 1, 1, 0),
    (1, 1, 1, 0),
    (0, 0, 0, 1),
    (1, 0, 0, 1),
    (0, 1, 0, 1),
    (1, 1, 0, 1),
    (0, 0, 1, 1),
    (1, 0, 1, 1),
    (0, 1, 1, 1),
    (1, 1, 1, 1),
]

initial_state_ls2 = [
    (0, 0, 0, 0),
    (1, 1, 1, 1),
]


def plotPhaseDiagram(initial_states, m1, m2, k1, k2, omega):
    t = np.linspace(0, 100, 1000)
    f = lambda x: np.sin(x*omega)
    for i, init in enumerate(initial_states):
        f_args = (m1, m2, k1, k2, f)
        sol = odeint(func, init, t, args = f_args)
        plt.plot(sol[:, 0], sol[:, 1])
        plt.grid()
        plt.show()


t, x1, x2 = computeODE((0, 0, 0, 0), 1, 1, 1, 1, 1)
plotODE(t, x1, x2)

initial_states = [
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1)
]


plotPhaseDiagram(initial_states, 1, 1, 1, 1, 0.5)
