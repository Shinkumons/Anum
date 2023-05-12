import scipy
import numpy
import math
from scipy import integrate
import matplotlib.pyplot as plt
import formule1
import sympy


def Cal(S0):
    """entrée : SO = condition initial de S
    Retourne l'erreur relative entre log(S0) et 
    (a / b)*(1 - S_infini) + log(S_infini)
    si log(S0) != 0,
    retourne l'erreur absolue sinon

    """
    Gauche = numpy.log(S0)
    S_inf = formule1.resolve(a, b, S0, 1-S0, [0, 1000])[1][0]
    Droite = (a / b)*(1 - S_inf) + numpy.log(S_inf)
    print("Le ln de S0 est :", Gauche)
    if Gauche ==0:
        return abs(Gauche-Droite)
    return abs((Droite - Gauche)/Gauche)

if __name__ == "__main__":
    a = 0.9
    b = 0.3
    S0 = 0.005
    while round(S0, 3) <= 1.0:
        erreur = Cal(round(S0, 3))
        print("Erreur  :", erreur,"|",  "S0 =", round(S0, 3))
        print("---------------------------------------------")
        S0 += 0.005
