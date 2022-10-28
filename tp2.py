import numpy as np
from numpy import linalg as la
import math
import matplotlib.pyplot as plt

t = 2
mat_eq = np.array([[math.sin(2*t), -1],
        [2*math.sin(t), math.cos(t)]])

mat_b = np.array([2*math.cos(t)*abs(math.sin(t))-1, 2*abs(math.sin(t)+math.cos(t))])

sol = la.solve(mat_eq, mat_b)

a, b = -1000, 1000
x, y = [], []
fig, axes = plt.subplots()

for i in np.linspace(a, b, 100000):
    mat_eq = np.array([[ math.sin( 2 * i ), -1],
        [2 * math.sin(i), math.cos( i )]])

    mat_b = np.array([2 * math.cos( i ) * abs( math.sin( i )) - 1, 2 * abs( math.sin( i )) + math.cos( i )])

    sol = la.solve(mat_eq, mat_b)

    x.append(sol[0])
    y.append(sol[1])

axes.grid()
axes.plot(x, y)
plt.show()
