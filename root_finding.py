import math
import unittest
import scipy.optimize
import matplotlib.pyplot as plt

# Q1 A
def is_good_enough(a, b, cnt, fan, fbn):
    """Return if a and b are almost equals"""

    return abs(a-b) <= 1e-15

def is_good_enough_newton(a, b, cnt, xn, f):
    """Return if a and b are almost equals (for newton)"""

    return abs(f(xn)) <= 1e-15

def bissection(f, a: float, b: float, max_iter=100, stop_cond=is_good_enough,
               is_plotting=False):
    """
    :return: root of f in the bounds of ]a, b[ using bissection method
    :param f: function
    :param a: lower bound for bissection
    :param b: upper bound for bissection
    :param max_iter: (optionnal) max number of iteration. If the iterations 
    exceed, will return None
    :param stop_cond: (optionnal) Function that take a, b, iteration number, 
    f(an) and f(bn) that stop the iteration as you want, currently return if 
    |a-b| <= 10^-15
    :param is_plotting: if True, return a tuple with the root and a list of the
    successives guesses of the bissection algorithm. otherise, return only the 
    root
    """
    if is_plotting:
        guess_root = []

    # Errors Management:  test if values are NaN or inf
    if math.isnan(a):
        raise ValueError("Parameter 'a' is NaN")
    if math.isnan(b):
        raise ValueError("Parameter 'b' is NaN")
    if math.isinf(a):
        raise ValueError("Parameter 'a' is inf")
    if math.isinf(b):
        raise ValueError("Parameter 'b' is inf")
    
    
    fan, fbn = f(a), f(b)
    if math.isnan(fan):
        raise ValueError("f is not defined in a")
    if math.isnan(fbn):
        raise ValueError("f is not defined in b")

    if f(a)*f(b) >= 0:
        raise ValueError(f"f(a)f(b) must be less than 0, current {f(a)*f(b)}")

    if fan == 0: return a
    if fbn == 0: return b

    
    # Permute a and b if needed
    if fan > 0 and fbn < 0:
        a, b, fan, fbn = b, a, fbn, fan
    xn = (a+b)/2
    # Loop algorithm
    cnt = 0
    while not stop_cond(a, b, cnt, fan, fbn) and cnt < max_iter:        
        xn = (a+b)/2
        if is_plotting:
            guess_root.append(xn)
        fxn = f(xn)
        if fxn > 0:
            b = xn
            fbn = fxn
        elif fxn < 0:
            a = xn
            fan = fxn
        else:
            if is_plotting:
                return (xn, guess_root[:-1])
            else:
                return xn
        cnt += 1
    if cnt >= max_iter:
        return None
    else:
        if is_plotting:
            return (xn, guess_root[:-1])
        else:
            return xn

def custom_stop_cond(a, b, cnt, fan, fbn):
    return False

class BissectUnitTesting(unittest.TestCase):
    def testBissect(self):
        self.assertAlmostEqual(bissection(lambda x: x**2 - 2, 0, 2),
                               math.sqrt(2))
        self.assertAlmostEqual(bissection(lambda x: x, -math.sqrt(2), 1), 0)
        self.assertIsNone(bissection(lambda x: -1 if x<0 else 1, -1, 1,
                                     stop_cond=custom_stop_cond))
        self.assertRaises(ValueError, bissection, lambda x: x**2, math.nan, 1)
        self.assertRaises(ValueError, bissection, lambda x: x**2, 1, math.nan)
        self.assertRaises(ValueError, bissection, lambda x: x**2, math.inf, 1)
        self.assertRaises(ValueError, bissection, lambda x: x**2, 1, math.inf)

def scipy_brentq(f, a, b):
    scipy.optimize.brentq(
        f,
        a,
        b,
        full_output=True,
        disp=True
    )

def embedded_plot(datas: list, root: float):
    fig = plt.figure()
    error_data = [abs(root-i) for  i in datas]
    axes = fig.add_axes([0.1, 0.1, 0.98, 0.98])
    y = list(range(len(datas)))
    axes.set_yscale('log')
    plt.plot(y, error_data)
    plt.grid()
    plt.show()


def embedded_plot_bissect_x_newton(data1: list, data2: list, root: float):
    fig = plt.figure()
    error_data1 = [abs(root-i) for  i in data1]
    error_data2 = [abs(root-i) for  i in data2]
    axes = fig.add_axes([0.1, 0.1, 0.98, 0.98])
    y1 = list(range(len(data1)))
    y2 = list(range(len(data2)))
    axes.set_yscale('log')
    plt.plot(y1, error_data1)
    plt.plot(y2, error_data2)
    plt.grid()
    plt.show()

    
def newton(f, f_der, a, b, stop_cond=is_good_enough_newton, max_iter=150, is_plotting=False):

    if is_plotting:
        guesses = []

    if math.isnan(a): raise ValueError("a parameter is NaN")
    if math.isnan(b): raise ValueError("b parameter is NaN")
    if math.isinf(a): raise ValueError("a parameter is infinite")
    if math.isinf(b): raise ValueError("n parameter is infinite")

    
    fa, fda = f(a), f_der(a)
    fb, fdb = f(b), f_der(b)

    if fdb == 0: raise ValueError("Derivative is null in b")
    if fda == 0: raise ValueError("Derivative is null in a")

    if math.isnan(fa): raise ValueError("f not defined in a")
    if math.isnan(fb): raise ValueError("f not defined in b")
    if math.isnan(fda): raise ValueError("f' not defined in a")
    if math.isnan(fdb): raise ValueError("f' not defined in b")

    if math.isinf(fa): raise ValueError("f is infinite in a")
    if math.isinf(fb): raise ValueError("f is infinite in b")
    if math.isinf(fda): raise ValueError("f' is infinite in a")
    if math.isinf(fdb): raise ValueError("f' is infinite in b")

    if fa*fb > 0: raise ValueError(f"No roots guarenteed between {a} and {b}, f(a)*f(b) > 0")
    
    xna, xnb = a-(fa/fda), b-(fb/fdb)
    
    if abs(f(xna)) <= abs(f(xnb)):
        xn = xna
    else:
        xn = xbn
    cnt = 0
    while not stop_cond(a, b, cnt, xn, f) or cnt >= max_iter:
        fxn, fdxn = f(xn), f_der(xn)
        xn -= fxn/fdxn
        if is_plotting:
            guesses.append(xn)
        if f(xn) == 0:
            if is_plotting:
                return xn, guesses
            else:
                return xn
        cnt += 1
    if cnt<max_iter:
        if is_plotting:
            return xn, guesses
        else:
            return xn
    else:
        return None
    
if __name__ == "__main__":
    print(scipy.optimize.brentq(lambda x: x**3 - 2, 1, 10))
    root1, guesses1 = newton(lambda x: x**3 -2, lambda x: 3*(x**2), 1, 10, is_plotting=True)
    
    root2, guesses2 =  bissection(lambda x: x**3 - 2, 1, 10, is_plotting=True)
    
    embedded_plot_bissect_x_newton(guesses1, guesses2, 2**(1/3))

    """
    On s'aperçois que la méthode de la bissection a une erreur qui décroit exponentiellement vers 0 (ce qui
    se traduit par une droite sur une échelle logarythmique) avec un graphe similaire a 1/2^n.
    
    On s'aperçois aussi que la méthode de Newton décroit exponentiellement vers 0 malgrès l'échelle logarythmique.
    Donc methode de Newton est beaucoup plus rapide que la bissecton.
    """
    
    unittest.main()
