import scipy
import numpy
import math
import sys
from scipy import integrate
import matplotlib.pyplot as plt
import formule1
import q8


def Precisionb2():
    """retourne le b minimal tq pour tout t, I(t)<=0.15
    """
    return scipy.optimize.brentq(fonction, 0, 1)

def maxI2(b):
    """entrée : b = Le pourcentage de personnes infectées qui deviennent
    immunisées
    retourne max I avec brentq
    """
    def fct(t):
        var = formule1.resolve(a, b, 1-I0, I0, [0, t])
        return deriveI(var, 1, b)
    var = formule1.resolve(a, b, 0.99, 0.01, [0, 50])
    if deriveI(var, 0, b) < 0 and deriveI(var, 1, b) < 0:
        return var[0][1]
    t = scipy.optimize.brentq(fct, 0, 50)
    var2 = formule1.resolve(a, b, 0.99, 0.01, [0, t])
    return var2[1][1]

def fonction(b):
    """ entrée : b = Le pourcentage de personnes infectées qui deviennent
    immunisées
    return maximum des I-0.15"""
    return maxI2(b) - 0.15


def deriveI(var, i, b):
    """entrée : var = tableau à n ligne 2 colonnes
                i   = indice ligne du tableau var
                b  = même chose que précédemment
    Retourne derive de I
    """
    return (a*var[i][0]- b) * var[i][1]

if __name__ == "__main__":
    t = numpy.linspace(0, 50, 100)
    I0 = 0.01
    S0 = 0.99
    a = 0.9
    b = 1
    print(Precisionb2())
