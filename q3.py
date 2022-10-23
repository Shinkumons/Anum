import scipy.optimize
import math

def is_zero(a):
    """This function return if a value is close enought to 0.
    for the q3 function.
    :param a: a float
    """

    return abs(a) <= 1e-10

def q3(a: float, y: float):
    """Prints pre-image of y by the function f(x) = x + a*sin(x) on the standard
    output
    :param a: a float, gives the range of sin(x).
    :param y: a float, the value of which you want the pre-images
    """

    # Errors managements
    if a<0: raise ValueError("a must be greater or equals to 0")
    if math.isNaN(a): raise ValueError("a is NaN")
    if math.isinf(a): raise ValueError("a is infinite")
    if math.isNaN(y): raise ValueError("y is NaN")
    if math.isinf(y): raise ValueError("y is infinite")

    # functions definitions
    f = lambda x: x + a*math.sin(x) - y
    df = lambda x: 1 + a*math.cos(x)

    # lower and upper bounds where we can find roots
    r_1 = y-a
    r_2 = y+a

    if a < 1:
        # management of the case where we cannot divide into intervals such that
        # the function is unimodal on these intervals

        # In particular, the function is strictly increasing (and thus
        # injective) because its derivative is > 0
        print(scipy.optimize.brentq(lambda x: x + a*math.sin(x) - y, r_1, r_2))

    else :
        # Otherwise we divide the search interval into intervals of length 2*pi
        # such that the function is unimodal on this interval.

        # find a particular root of the derivative (we need that a >= 1
        # otherwise, the derivative don't have roots)
        x_0 = -math.acos(-1/a)
        # calculate lowest and greater k value for x_0 + 2*k*pi
        k_1 = math.floor((r_1 - x_0) / (2*math.pi))
        k_2 = math.ceil((r_2 - x_0)/(2*math.pi))


        # Root searching in the intervals
        for i in range(k_1, k_2):
            # bounds of the intervals
            b_1, b_2 = x_0 + i*2*math.pi, x_0 + (i+1)*2*math.pi

            m = -x_0 + 2*i*math.pi
            fb_1, fm, fb_2 = f(b_1), f(m), f(b_2)

            if fb_1 < 0:
                print(scipy.optimize.brentq(f, b_1, m))

            elif is_zero(fb_1):
                print(b_1)

            if fb_2 < 0:
                print(scipy.optimize.brentq(f, m, b_2))
