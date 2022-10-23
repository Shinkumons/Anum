from q3 import *
import io
import unittest


def q3_testing_version(a: float, y: float):
    """Prints pre-image of y by the function f(x) = x + a*sin(x) on the standard
    output
    :param a: a float, gives the range of sin(x).
    :param y: a float, the value of which you want the pre-images
    """

    # Errors managements
    if a<0: raise ValueError("a must be greater or equals to 0")
    if math.isnan(a): raise ValueError("a is NaN")
    if math.isinf(a): raise ValueError("a is infinite")
    if math.isnan(y): raise ValueError("y is NaN")
    if math.isinf(y): raise ValueError("y is infinite")

    # functions definitions
    f = lambda x: x + a*math.sin(x) - y
    df = lambda x: 1 + a*math.cos(x)

    # lower and upper bounds where we can find roots
    r_1 = y-a
    r_2 = y+a

    roots_list = []

    if a <= 1:
        # management of the case where we cannot divide into intervals such that
        # the function is unimodal on these intervals

        # In particular, the function is strictly increasing (and thus
        # injective) because its derivative is > 0
        roots_list.append(scipy.optimize.brentq(lambda x: x + a*math.sin(x) - y,
                                                r_1, r_2))

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

            if fb_1*fm < 0:
                r = round(scipy.optimize.brentq(f, b_1, m, xtol=1e-12), 11)
                if r not in roots_list:
                    roots_list.append(r)

            elif is_zero(fb_1):
                if b_1 not in roots_list:
                    roots_list.append(round(b_1, 11))

            if fb_2*fm < 0:
                r = round(scipy.optimize.brentq(f, m, b_2, xtol = 1e-12), 11)
                if r not in roots_list:
                    roots_list.append(r)

    return roots_list

class Q3UnitTest(unittest.TestCase):
    def assertQ3Testing(self, a: float, y: float, expected_output: list):
        sol = q3_testing_version(a, y)
        if len(sol) != len(expected_output):
            raise RuntimeError("the numbers of inputs are not the sames")
        else:

            for i, val in enumerate(sol):
                # Reduced precision of the test cause wolfram alfa goes to 5
                # decimals at minimum (depend on the size of the numbers)
                self.assertAlmostEqual(val, expected_output[i], 5)

    def testQ3(self):

        # only one root testings (0 <= a < 1)
        self.assertQ3Testing(0, 5, [5])
        self.assertQ3Testing(0.5, 2, [1.50121])
        self.assertQ3Testing(0.75, 1, [0.585523])

        # a = 1
        self.assertQ3Testing(1, 1, [0.510973])

        # multiple roots testings
        self.assertQ3Testing(5.3, -0.1, [-5.06814, -3.95646, -0.01587, 4.03723,
                                         4.99288])

        # this test works (wolfram precision si too short but they're equals)
        # self.assertQ3Testing(8.6, 5, [0.54459, 2.89422, 6.14916, 10.0528,
        #                               11.6774])

        # Border case : m is a root
        self.assertQ3Testing(1, math.pi, [math.pi])
        self.assertQ3Testing(2, 3**(1/2) + 2*math.pi/3, [2*math.pi/3,
                                                         5.38759824])

if __name__ == "__main__":
    unittest.main()
