from math import *
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg
from scipy.optimize import *

def remez_cond(p_n, p_n1, a, b):
    p_diff = lambda x: p_n(x) - p_n1(x)
    eval_p_diff = p_diff(np.linspace(a, b, ceil(abs(100*(a-b)))))

    p_diff_m = abs(max(eval_p_diff))
    return p_diff_m <= 1e-5*(abs(a)+abs(b)) + 1e-10

def min_max_search(err, roots, ifm):
    min_max = []
    for i in range(len(roots) - 1):
        x, y = roots[i], roots[i+1]
        if ifm:
            maxi, val = 0, 0
            ls = np.linspace(x, y, ceil(abs(100*(y-x))))
            for i, value in enumerate(ls):
                if err(value) > val:
                    maxi, val = i, err(value)
            min_max.append(ls[maxi])
        else:
            mini, val = 0, 0
            ls = np.linspace(x, y, ceil(abs(100*(y-x))))
            for i, value in enumerate(ls):
                if err(value) < val:
                    mini, val = i, err(value)
            min_max.append(ls[mini])

        ifm = not ifm

    return min_max
            

def remez(f, n, a, b):

    cheb_node = [ 0.5*(a+b) + 0.5*(b-a) * cos((((2*k)-1)/(2*(n+1)))*pi) for k in range(1,n+2)][::-1]

    p_0 = interpolate.lagrange(np.array(cheb_node), np.array([f(x) for x in cheb_node]))# , kind='quadratic')
    err_0 = lambda x: f(x) - p_0(x)
    
    tem_point = (cheb_node[0] + cheb_node[1])/2
    if err_0(tem_point) < 0 :
        is_first_max = False
    else:
        is_first_max = True

    x = np.linspace(a, b, 1000)
    y = err_0(x)
        
    x_0 = np.array([a] + min_max_search(err_0, cheb_node, is_first_max) + [b])
    
    lin_sys = [[x_0[i]**j for j in range(n+1)] + [(-1)**i] for i in range(len(x_0))]
    
    sol = linalg.solve(np.array(lin_sys), f(x_0))

    p_1 = np.poly1d(sol[::-1][1:])

    p_n = p_0
    p_n1 = p_1
    e_n = err_0
    
    while not remez_cond(p_n, p_n1, a, b):
        e_n1 = lambda x: f(x) - p_n1(x)
        
        roots = [brentq(e_n1, x_0[i], x_0[i+1]) for i in range(len(x_0)-1)]

        tem_point = (roots[0] + roots[1])/2
        if err_0(tem_point) < 0 :
            is_first_max = False
        else:
            is_first_max = True

        minmax = min_max_search(e_n1, roots, is_first_max)
        x_0 = np.array([a] + minmax + [b])
        lin_sys = [[x_0[i]**j for j in range(n+1)] + [(-1)**i] for i in range(len(x_0))]
        sol = linalg.solve(np.array(lin_sys), f(x_0))
        p_n = p_n1
        p_n1 = np.poly1d(sol[::-1][1:])

    return x_0, p_n1

f = lambda x: e**x
x_i, p_i = remez(f, 2, -2, 2)
x = np.linspace(-2, 2, 1000)
plt.plot(x, f(x), color="tab:blue")
plt.plot(x, p_i(x), color="tab:red")
plt.grid()
plt.show()

