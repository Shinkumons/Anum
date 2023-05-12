import scipy
import numpy
import math
import sys
from scipy import integrate
import matplotlib.pyplot as plt
import random
import sympy
"""
Le programme prend un argument lors de son execution : 
python formule1.py args
args est un naturel.

a = le nombre moyen d’individus qui ont un contact suffisant
pour transmettre l’infection par
unité de temps ;
b = Le pourcentage de personnes infectées qui deviennent immunisées
par unité de temps
b represente alpha dans nos équations différentielles.
S0 = condition initiale de S ( taux de la population saine )
I0 = condition initiale de I ( taux de la population infectée )
t = un intervale de temps

Ces paramètres sont utilisé comme entrée dans les différentes méthodes
ci-dessous.
"""

def solve(a, b, S0, I0, t):
    """ entrée : voir ci-dessus
    Retourne la fonction solution de l'équation différentielle
    """
    def fonction(S0, I0, t):
        return integrate.odeint(mush, [S0, I0], t)
    def mush(y, t):
        S, I = y
        return [-a*S*I, (a*S - b)*I]
    return fonction

def resolve(a, b, S0, I0, t):
    """ Retourne la fonction solution en t
    """
    var = solve(a, b, S0, I0, t)
    return var(S0, I0, t)

def converter(Tab):
    """ Entrée : tableau à n ligne 2 colonnes
    sortie : tableau 2 lignes n colonnes
    tableau à n ligne 2 colonnes qu'on transforme en un tableau
    à 2 ligne et n colonnes
   """
    S = []
    I = []
    for i in range(len(Tab)):
        S += [Tab[i][0]]
        I += [Tab[i][1]]
    return [S, I]

def Plsvar(n):
    """entrée : naturel n
    Permet de créer n solutions de l'équation différentielle avec
    des conditions initiales créées de manière aléatoire
    sortie : dictionnaire possédant les solutions des différentes
    équations.
    """
    dico = {}
    for i in range(n):
        S0 = float("%.3f"%random.random())
        I0 = float("%.3f"%random.random())
        while S0 + I0 >1:
            S0 = float("%.3f"%random.random())
            I0 = float("%.3f"%random.random())
        dico[i] = resolve(a, b, S0, I0, t)
    return dico



if __name__ == "__main__":
    input = int(sys.argv[1]) #nombre de conditions initiales
    a = 0.9
    b = 0.3
    t = numpy.linspace(0, 50, 1000)
    var = Plsvar(input)
    plt.figure(figsize =(12, 9))
    for i in range(input):
        var[i] = converter(var[i])
        plt.plot(var[i][0], var[i][1], label = ('S0 : '+ str(var[i][0][0]),
        'I0 : ' + str(var[i][1][0])))
    plt.xlabel("S")
    plt.ylabel("I")
    plt.plot([0, 1], [1, 0], color = "black")
    plt.plot([0, 1], [0, 0], color = "black")
    plt.plot([0, 0], [1, 0], color = "black")
    plt.legend(bbox_to_anchor=(1.10, 0.5), loc = "center right")
    plt.show()
