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

def plot_array_ODE(initial_states, k1, k2, m1, m2, omega):
    for init in initial_states:
        t, x1, x2 = computeODE(init, k1, k2, m1, m2, omega)
        plotODE(t, x1, x2)

initial_states = [
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (0, 0, 0, 0.1)
]
omega = 1
plotDiagrams(initial_states, 1, 1, 1, 1, omega)
