import scipy
import numpy
import math
import sys
from scipy import integrate
import matplotlib.pyplot as plt
"""
   g représente gamma
   n représente eta
   d représente delta
   dans les formules (3a)-(3c)
   a, b, S0, I0, t sont identique à ceux de formule1
   g = fraction de la population infecté qui recoit le traitement
   par unité de temps 
   n = taux de guérison dû au traitement
   d = taux de reduction de l'infectiosité d'une personne sous traitement
   T0 = condition initiale de T (taux de la population sous traitement)
   
   ces paramètres sont utilisées comme entrées dans les méthodes ci-dessous
"""

def solve(a, b, g, n, d, S0, I0, T0, t):
    """Retourne la fonction solution de l'équation différentielle
    """
    def fonction(S0, I0, T0, t):
        return integrate.odeint(mush, [S0, I0, T0], t)
    def mush(y, t):
        S, I, T = y
        return [-a*S*(I+d*T), a*S*(I+d*T) - (b+g)*I, g*I-n*T]
    return fonction

def resolve(a, b, g, n, d, S0, I0, T0, t):
    """Retourne la fonction solution en t
    """
    var = solve(a, b, g, n, d, S0, I0, T0, t)
    return var(S0, I0, T0, t)

if __name__=="__main__":
    t = numpy.linspace(0, 50, 100)
    a = 0.9
    b = 0.1
    d = 0.3
    n = 0.5
    g = 0.2
    Sol = resolve(a, b, g, n, d, 0.99, 0.01, 0, t)
    plt.figure(figsize =(12, 8))#taille
    plt.plot(t, Sol, label = "1")
    plt.show()
