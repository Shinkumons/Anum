import formule1
import scipy
import numpy
import math
import sys
from scipy import integrate
import matplotlib.pyplot as plt


def maxI2(I0):
    """entrée : I_0 = condition initiale de I
    Retourne max I avec brentq
    """
    def fct(t):
        var = formule1.resolve(a, b, 1-I0, I0, [0, t])
        return deriveI(var, 1)
    var = formule1.resolve(a, b, 1-I0, I0, [0, 50])
    if (a*var[0][0]-b)*var[0][1] < 0 and (a*var[1][0]-b)*var[1][1] < 0:
        return var[0][1]
    t = scipy.optimize.brentq(fct, 0, 50)
    var2 = formule1.resolve(a, b, 1-I0, I0, [0, t])
    return var2[1][1]


def Precision2():
    """retourne le I0 minimal tq pour tout t, I(t) <=0.4
    """
    return scipy.optimize.brentq(fonction, 0, 0.4)

def fonction(I0):
    """ entrée : I_0
    Retourne le maximum de I-04
    """
    return maxI2(I0) - 0.4


def deriveI(var, i):
    """entrée : var = tableau à n ligne 2 colonnes
                i   = indice ligne du tableau var
    Retourne la dérive de I en un certain temps que i représente
    """
    return (a*var[i][0]- b) * var[i][1]

if __name__=="__main__":
    t = numpy.linspace(0, 50, 100000)
    a = 0.9
    b = 0.3
    Imin = Precision2()
    print(Precision2())
    Sol = formule1.resolve(a, b, 1-Imin, Imin, t)
    plt.figure()
    plt.plot(t, Sol, label = "I0")
    plt.show()