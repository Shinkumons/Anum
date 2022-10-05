import math
import unittest
import scipy.optimize
import matplotlib.pyplot as plt

# Q1 A
def is_good_enough(a, b, cnt, fan, fbn):
    """Return if a and b are almost equals"""

    return abs(a-b) <= 1e-15

def bissection(f, a, b, max_iter=100, stop_cond=is_good_enough,
               is_plotting=False):
    """
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
    
    Return a root of f in the bounds of ]a, b[ using bissection method
    """
    if is_plotting:
        guess_root = []

    # Errors Management:  test if values are NaN or inf
    if f(a)*f(b) >= 0:
        raise ValueError(f"f(a)f(b) must be less than 0, current {f(a)*f(b)}")
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
                return (xn, guess_root)
            else:
                return xn
        cnt += 1
    if cnt >= max_iter:
        return None
    else:
        if is_plotting:
            return (xn, guess_root)
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
        
if __name__ == "__main__":
    print(scipy.optimize.brentq(lambda x: x**3 - 2, 0, 50))
    root, guesses =  bissection(lambda x: x**3 - 2, 0, 2, is_plotting=True)
    guesses = [abs(root-i) for  i in guesses]
    fig = plt.figure()
    axes = fig.add_axes([0.1, 0.1, 0.98, 0.98])
    y = list(range(len(guesses)))
    axes.set_yscale('log')
    plt.plot(y, guesses)
    plt.grid()
    plt.show()
    unittest.main()
