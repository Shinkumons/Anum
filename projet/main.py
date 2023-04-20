"""Projet épidémiologie

Usage:
  projet_epi.py -c <a> <alpha>
  projet_epi.py -h
  projet_epi.py -max <a> <alpha> <i>

Options:
  -h            montre cette fenêtre
  -c            calcule l'EDO en fonction de <a> et <alpha> et affiche le triangle des solutions
  -max           calcule le nombre d'infecté minimal pour atteindre un seuil d'infection i
"""

from scipy.integrate import odeint
from scipy.optimize import brentq
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from docopt import docopt

# fonction représenant l'EDO
# dS/dx = -a*S*I
# dI/dx = (a*S-alpha)*I
def equations(y, t, a, alpha):
    '''
    Calcule la dérivée de S et I en fonction de y et t des équations de paramètres a et alpha.

    :param y: un couple de float représentant les solutions précédentes des équations.
    :param t: un nombre représentant l'évaluation en t de S et I.
    :param a: le coefficient représentant le taux d'infection du système, doit être > 0.
    :param alpha: le coefficient représentant le taux de rémission du système, doit être > 0.
    :return: un couple représentant la dérivée de S et I au moment t. 
    '''
    s, i = y
    equations = [
        -a*s*i,
        (a*s - alpha)*i
    ]
    return equations

# Calcul de l'EDO en fonction des conditions initiales s0, i0 et des paramètre a et alpha
def computeODE(initial_state, a, alpha):
    '''
    Calcule une solution du modèle SIR en fonction des conditions initiales s_0 et i_0 en t=0 et les coefficients a et alpha.
    
    :param initial_state: tuple s_0, i_0 représentant les conditions initiales des équations différentielles.
    :param val_a: valeur de a, un nombre a > 0.
    :param val_alpha: valeur de alpha, un nombre alpha > 0.
    :return: 3 listes, la liste des temps t, la listes des images de t par S et une liste des images de t par I. 
    '''
    
    t = np.linspace(0, 40, 1000)
    f_args = (a, alpha)
    sol = odeint(equations, initial_state, t, args = f_args)
    return (t, sol[:, 0], sol[:, 1])
  
    
def findMinimumInfected(a, alpha, find):
    f = lambda x: 1 + (alpha/a)*(-1 + np.ln(alpha/(a*(1-x)))) - find

    sol = brentq(f, 0, 1)
    return sol

def calculateTriangle(val_a, val_alpha):
    """
    Calcule et affiche une multitude de solution en fonction de a et alpha.

    :param val_a: valeur de a, un nombre a > 0.
    :param val_alpha: valeur de alpha, un nombre alpha > 0.
    :return: 
    """
    

    # Représentation du triangle pour la fonction plt.fill() de matplotlib.
    triangle_x = [0, 1, 0]
    triangle_y = [0, 0, 1]

    # Plot du triangle T.
    plt.fill(triangle_x, triangle_y, alpha = 0.1)
    
    # Borne inférieure du nombre d'infecté initial, le cas ou le nombres d'infecté étant constant, il n'est pas intéressant
    # donc on la met proche de 0.
    init_inf = 0.01

    # Nombres de subdivisions, /!\ augmenter ce paramètre augmente le temps de calcul.
    mesh_subdiv = 20
    # Liste des conditions initiales au bord du triangle donnant le graphe le plus satifaisant et clair
    init_cond = [(i, init_inf) for i in np.linspace(0.5, 1-init_inf, mesh_subdiv//2)]
    init_cond += [(1-i, i) for i in np.linspace(0, 1, mesh_subdiv)]
    
    for i, j in init_cond:
        time, x_sol, y_sol = computeODE((i, j), val_a, val_alpha)
        plt.plot(x_sol, y_sol, color="red")
        
        # Plot des flèches indiquant le sens des solutions (pas ideal mais faute d'avoir mieux pour l'instant).
        for arr_p in range(1, 4):
            k = math.floor(len(x_sol)* arr_p/10)
            arrow_vect = (x_sol[k+1],
                          y_sol[k+1],
                          x_sol[k+1] - x_sol[k],
                          y_sol[k+1] - y_sol[k])
            plt.arrow(*arrow_vect, shape='full', lw=0,
                      length_includes_head=True, head_width=0.01, color="red")

    plt.grid()
    plt.show()


# TODO remplacer ce truc ignoble par la documentation docopt !
if __name__ == "__main__":
    arguments = docopt(__doc__, version='0.1.1rc')
    print(arguments)

    if arguments['-c']:
        try:
            a, alpha = float(arguments['<a>']), float(arguments['<alpha>'])
            calculateTriangle(a, alpha)
        except:
            print("a et alpha doivent être des valeurs réelles")

    if arguments['-max']:
        try:
            a, alpha, i = float(arguments['<a>']), float(arguments['<alpha>']), float(arguments['<i>'])
            print(findMinimumInfected(a, alpha, i))
        except:
            print("a et alpha et i doivent être des valeurs réelles")
