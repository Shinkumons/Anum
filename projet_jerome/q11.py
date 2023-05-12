import scipy
import numpy
import math
import sys
from scipy import integrate
import matplotlib.pyplot as plt
import formule2


def maxI2(g):
    """entrée : g = fraction de la population infecté qui recoit
   le traitement par unité de temps 
    Retourne max de I avec brentq
    """
    def fct(t):
        var = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, t])
        return deriveI(var, 1, g)
    var = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, 50])
    if deriveI(var, 0, g) < 0 and deriveI(var, 1, g) < 0:
        return var[0][1]
    t = scipy.optimize.brentq(fct, 0, 50)
    var2 = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, t])
    print(var2)
    return var2[1][1]

def Precisiong2():
    """retourne le gamma minimal tq pour tout t, I(t)<=0.15
    """
    return scipy.optimize.brentq(fonction, 0, 1)

def fonction(g):
    """entrée : g = fraction de la population infecté qui recoit
   le traitement par unité de temps 
    Retourne le maximum de I - 0.15
    """
    return maxI2(g) - 0.15

def deriveI(var, i, g):
    """entrée : var : tableau à n lignes et 3 colonnes
                i : indice de la ligne
                g : même qu'au dessus
    Derive de I"""
    return (a*var[i][0])*(var[i][1] + d*var[i][2]) - (b+g)*var[i][1]


#-------------------------

def deriveT(var, i, g):
    """entrée : var : tableau à n lignes et 3 colonnes
                i : indice de la ligne
                g : même qu'au dessus
    derive de T"""
    return g*var[i][1] - n*var[i][2]


def maxT2(g):
    """Retourne max de T avec brentq
    """
    def fct(t):
        var = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, t])
        return deriveT(var, 1, g)
    var = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, 50])
    if deriveT(var, 0, g) < 0 and deriveT(var, 1, g) < 0:
        return var[0][2]
    t = scipy.optimize.brentq(fct, 0, 50)
    var2 = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01,T0, [0, t])
    return var2[1][2]


if __name__=="__main__":
    t1 = numpy.linspace(0, 50, 100)
    a = 0.9
    b = 0.1
    d = 0.3
    n = 0.5
    I0 = 0.01
    S0 = 0.99
    T0 = 0
    g = Precisiong2()
    tab = []
    tab2 = []
    print("valeur maximum de gamma :" ,g)
    print("maximum de T avec ce gamma :", maxT2(g))
    #Sol = formule2.resolve(0.9, 0.1, g, 0.5, 0.3, 0.99, 0.01, T0, t1)
    #plt.figure()
    #plt.plot(t1, Sol, label = "1")
    #plt.plot(tab, color = "black")
    #plt.plot(tab2, color = "red")
    #plt.show()
