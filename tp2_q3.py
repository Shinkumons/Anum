import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def indic(x):
    if 0 <= x < math.pi:
        return 1
    else:
        return 0

def indic2(x):
    if 0<= x < math.pi:
        return 0.5
    else:
        return -0.5

def f_test(x):
    return (x%math.pi)
    
def approx(f, n):
    Bk = [42]
    Ck = []

    actual_ck = (1/(math.pi*2)) * integrate.quad(lambda x: f(x), 0, 2*math.pi, limit = 200)[0]
    Ck.append(actual_ck)
    
    for k in range(1, n+1):
        actual_bk = (1/math.pi) * integrate.quad(lambda x: np.sin(k*x)*f(x), 0, 2*math.pi, limit = 200)[0]
        Bk.append(actual_bk)
        actual_ck = (1/math.pi) * integrate.quad(lambda x: np.cos(k*x)*f(x), 0, 2*math.pi, limit = 200)[0]
        Ck.append(actual_ck)

    approx_f = lambda x: sum([Bk[i]*np.sin(i*x) for i in range(len(Bk))]) + sum([Ck[i]*np.cos(i*x) for i in range(len(Ck))])

    return approx_f, Bk, Ck

#approx_f, Bk, Ck = approx(indic, 200)
f = lambda x: (1/8)*x**3 - x**2 +2*x - 1
approx_f, Bk, Ck = approx(indic, 100)

x = np.linspace(0, 2*math.pi, 1000)
y1 = approx_f(x)

y2 = list(map(indic, x))

plt.plot(x, y1)
plt.plot(x, y2, '--')
plt.grid()
plt.show()

approx_f, Bk, Ck = approx(f, 100)

y1 = approx_f(x)
y2 = f(x)

plt.plot(x, y1)
plt.plot(x, y2, '--')
plt.grid()
plt.show()
