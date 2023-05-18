"""Projet épidémiologie
Usage:
  projet_epi.py -h
  projet_epi.py --q6 <a> <alpha>
  projet_epi.py --q8 <a> <alpha> <i>
  projet_epi.py --q9 <a> <alpha>
  projet_epi.py --q9 <a> <alpha> <nbr_test>
  projet_epi.py --q10 <a> <i0> <find>
  projet_epi.py --q11

Options:
  -h            montre cette fenêtre
  --q6          calcule l'EDO en fonction de <a> et <alpha> et affiche le
                triangle des solutions
  --q8          calcule le nombre d'infecté minimal pour atteindre un seuil
                d'infection i
  --q10
  --q11
"""
from scipy.integrate import odeint
from scipy.optimize import brentq
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from docopt import docopt

"""
Ce projet vise a étudier le comportement d'une épidémie au moyen
du modèle SIR (Susceptible Infectious Removed).

Nous prendrons donc des conventions de nommages pour la suite :

a := Le nombre moyen de personnes ayant un contact suffisant pour
transmettre l'infection par unités de temps
alpha := Le pourcentage de personnes infectée qui deviennent
immunisées par unité de temps
"""







# fonction représenant l'EDO
# dS/dx = -a*S*I
# dI/dx = (a*S-alpha)*I
def sir(y, t, a, alpha):
    '''
    Calcule la dérivée de S et I en fonction de y et t des équations de
    paramètres a et alpha.

    :param y: un couple de float représentant les solutions précédentes
              des équations.
    :param t: un nombre représentant l'évaluation en t de S et I.
    :param a: le coefficient représentant le taux d'infection du système,
              doit être > 0.
    :param alpha: le coefficient représentant le taux de rémission du système,
              doit être > 0.
    :return: un couple représentant la dérivée de S et I au moment t.
    '''
    s, i = y
    si = [
        -a*s*i,
        (a*s - alpha)*i
    ]
    return si

# Calcul de l'EDO en fonction des conditions initiales s0, i0 et des paramètre
# a et alpha
def computeODE(initial_state, a, alpha):
    '''
    Calcule une solution du modèle SIR en fonction des conditions initiales
    s_0 et i_0 en t=0 et les coefficients a et alpha.

    :param initial_state: tuple s_0, i_0 représentant les conditions initiales
                          des équations différentielles.
    :param a: valeur de a, un nombre a > 0.
    :param alpha: valeur de alpha, un nombre alpha > 0.
    :return: 3 listes, la liste des temps t, la listes des images de t par S et
             une liste des images de t par I.
    '''

    t = np.linspace(0, 40, 1000)
    f_args = (a, alpha)
    sol = odeint(sir, initial_state, t, args = f_args)
    return (t, sol[:, 0], sol[:, 1])


#Q8
def findMinimumInfected(a, alpha, find):
    """
    Calcule le taux minimum d'infecté initial tel que l'infection dépasse un
    seuil donné

    :param a: valeur de a, un nombre a > 0.
    :param alpha: valeur de alpha, un nombre alpha > 0.
    :param find: le seuil d'infection souhaité
    :return: la valeur minimale du taux d'infectés initial tel que l'infection
             dépasse le seuil
    """

    """
    Fonction représentant le pourcentage d'infectés maximums
    en fonction du pourcentage d'infectés initial et des facteurs
    a et alpha (voir évellopement mathématique dans le rapport)
    ici on cherche la valeur find et donc on la soustrait a la
    fonction"""
    f = lambda x: 1 + (alpha/a)*(-1 + np.log(alpha/(a*(1-x)))) - find

    # recherche de la racine
    sol = brentq(f, 0, find)
    return sol

def relative_error(a, b, x_tol, r_tol):
    return abs(a-b) <= r_tol * (abs(a)+abs(b)) + x_tol

def q9(a, alpha, nbr_test =  100):
    """Teste experimentalement la relation suivante

    \log(S_0) = \frac{\alpha}{a}(1-S_{\infty})\log(S_{\infty})

    :param a: valeur de a, un nombre a > 0.
    :param alpha: valeur de alpha, un nombre alpha > 0.
    :param nbr_test: (optionnel) nombre de tests effectués (valeur de base
                     fixée à 100)
    :return: Le résultat du test sous forme d'un pourcentage indiquant la
             proportion de test réussis
    """
    max_err = 0
    passed_test = 0
    # boucle sur plusieurs valeurs d'infectés initials
    for i in range(1, nbr_test+1):
        s0 = i/(nbr_test+1)
        i0 = 1-s0
        initial_state = (s0, i0)
        f_args = (a, alpha)
        # on estime que la valeur de 100_000_000_000 est
        # suffisante pour estimes S_inf
        t = np.linspace(0, 100_000_000_000, 100000)
        # résolution de l'équation différentielle
        sol = odeint(sir, initial_state, t, args = f_args)
        s = sol[:, 0]
        s_inf = s[-1]
        value = (a/alpha)*(1-s_inf) + np.log(s_inf)
        # test de l'erreur
        err = abs(value - np.log(s0))
        if err > max_err : max_err = err
        if relative_error(np.log(s0), value, 1e-5, 1e-8):
            passed_test += 1
    print(max_err)
    return passed_test/nbr_test

#Q10
def q10(a, i0, find):
    """Calcule le alpha minimal qu'il faudrais atteindre pour garantir que
    l'infection ne dépasse pas le seuil find en fonction du facteur a et i0

    :param a: valeur de a, un nombre a>0
    :param i0: le porcentage d'infectés initial
    :param find: le pourcentage d'infectés maximal souhaités
    :return: retourne le facteur alpha à atteindre
    """

    """Cette fonction indique le nombre d'infectés initial avec a et find fixé,
    find le pourcentage maximal d'infection souhaité (développement
    mathématique dans le rapport). On effectue ensuite une recherche
    de racine sur cette fonction avec brentq
    """
    f = lambda x: -(x/a)*math.pow(math.e, -(a/x)*(find-1)-1) + 1 - i0

    upper_bound = 1.5 * (a * (1-find))
    lower_bound = 0.5 * (a * (1-find))

    sol = brentq(f, lower_bound, upper_bound)
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

    # Borne inférieure du nombre d'infecté initial, le cas ou le nombres
    # d'infecté étant constant, il n'est pas intéressant donc on la met
    # proche de 0.
    init_inf = 0.01

    # Nombres de subdivisions, /!\ augmenter ce paramètre augmente le temps de
    # calcul.
    mesh_subdiv = 20
    # Liste des conditions initiales au bord du triangle donnant le graphe le
    # plus satifaisant et clair
    treshold = mesh_subdiv//2
    init_cond = [(i, init_inf) for i in np.linspace(0.5, 1-init_inf, treshold)]
    init_cond += [(1-i, i) for i in np.linspace(0, 1, mesh_subdiv)]

    for i, j in init_cond:
        time, x_sol, y_sol = computeODE((i, j), val_a, val_alpha)
        plt.plot(x_sol, y_sol, color="red")

        # Plot des flèches indiquant le sens des solutions (pas ideal mais faute
        # d'avoir mieux pour l'instant).
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


def solveSITR(a, alpha, gamma, eta, delta, initial_state, t):
    """
    Retourne une lambda fonction résolvant le système SITR pour les conditions
    initiales pour des paramètres a, alpha, gamma, eta et delta.

    :param a: taux d'infection
    :param alpha: taux de guérison naturelle
    :param gamma: taux de traitement
    :param eta: taux de guérison
    :param delta: taux de réduction de transmission
    :param initial_state: conditions initiales de l'edo
    :param t: itérable représentant le temps
    :return: une fonction lambda représentant l'edo de paramètres donnés dont
    les variables sont les conditions initiales
    """
    f_args = (a, alpha, gamma, eta, delta)
    def sitr(y, t, a, alpha, gamma, eta, delta):
        S, I, T = y
        sit = [
            -a*S*(I + delta*T),
            a*S*(I + delta*T) - (alpha+gamma)*I,
            gamma*I - eta*T
        ]
        return sit

    f=lambda initial_state, t: odeint(sitr,
                                   initial_state,
                                   t,
                                   args = f_args)
    return f


def evalSITR(a, alpha, gamma, eta, delta, initial_state, t):
    """
    Evaluation de la fonction retournée par solveSITR

    :param a: taux d'infection
    :param alpha: taux de guérison naturelle
    :param gamma: taux de traitement prodigué
    :param eta: taux de guérison
    :param delta: taux de réduction de transmission
    :param initial_state: conditions initiales de l'edo
    :param t: itérable représentant le temps
    :return: evaluation de la fonction retournée par solveSITR
    """
    ev = solveSITR(a, alpha, gamma, eta, delta, initial_state, t)
    return ev(initial_state, t)


def computeGamma(a, alpha, eta, delta, initial_state, obj):
    """
    Retourne le gamma minimum à atteindre pour que l'infection ne dépasse
    pas la valeur obj

    :param a: taux d'infection
    :param alpha: taux de guérison naturelle
    :param eta: taux de guérison
    :param delta: taux de réduction de transmission
    :param initial_state: conditions initiales de l'edo
    :param obj: objectif du taux d'infection maximal de l'épidémie
    :return: le gamma tel que le taux d'infecté maximal de l'épidémie
             est obj
    """
    def maxI(gamma):
        """
        fonction retournant la valeur maximale de I en fonction de gamma
        """
        def deriv(ev, i):
            # dérivée de I
            return a*ev[i][0]*(ev[i][1]+delta*ev[i][2]) - (alpha+gamma)*ev[i][1]

        def fct(t):
            # fonction retournant la dérivée de I(t)
            ev = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, t])
            return deriv(ev, 1)

        # On évalue aux bornes, si chacune est négative alors il n'y a pas de
        # solution
        ev = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, 40])
        if deriv(ev, 0) < 0 and deriv(ev, 1) < 0:
            return ev[0][1]
        # recherche d'une racine
        t = brentq(fct, 0, 40)
        # evaluation de I en son maximum
        ev2 = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, t])
        return ev2[1][1]
    # recherche de la racine de la fonction maxI
    to_max = lambda g: maxI(g)-obj

    return brentq(to_max, 0, 40)

def computeTMax(a, alpha, gamma, eta, delta, initial_state):
    """
    Retourne la valeur maximum de la fonction T

    :param a: taux d'infection
    :param alpha: taux de guérison naturelle
    :param gamma: taux de traitement prodigué
    :param eta: taux de guérison
    :param delta: taux de réduction de transmission
    :param initial_state: conditions initiales de l'edo
    :return: la valeur maximale de la fonction T de l'edo
    """

    
    # dérivée de T, issue de des equations (3a)-(3c)
    deriv = lambda ev, i: gamma*ev[i][1]-eta*ev[i][2]

    # fonction retournant la dérivée de T en evaluant l'EDO
    def fct_to_max(t):
        ev = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, t])
        return deriv(ev, 1)

    # On evalue aux deux bornes de recherches, si les deux sont négatives,
    # cela veut dire qu'il n'y a pas de racines
    ev = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, 40])
    if deriv(ev, 0) < 0 and deriv(ev, 1) < 0:
        return ev[0][2]
    # recherche de la racine de la dérivée pour trouver le maximum
    t = brentq(fct_to_max, 0, 40)
    # evaluation de l'edo en t
    ev2 = evalSITR(a, alpha, gamma, eta, delta, initial_state, [0, t])
    # retour de la fonction T(t) en t
    return ev2[1][2]



# Cette partie de code gère les arguments passés en ligne de commande pour plus
# de simplicité mais les fonctions de ce projets peuvent évidemment être
# importée dans un autre code
if __name__ == "__main__":
    arguments = docopt(__doc__, version='0.1.1rc')

    if arguments['--q6']:
        try:
            a, alpha = float(arguments['<a>']), float(arguments['<alpha>'])
            calculateTriangle(a, alpha)
        except:
            print("a et alpha doivent être des valeurs réelles")

    if arguments['--q8']:
        try:
            a, alpha = float(arguments['<a>']), float(arguments['<alpha>'])
            i = float(arguments['<i>'])
            print(findMinimumInfected(a, alpha, i))
        except:
            print("a et alpha et i doivent être des valeurs réelles")

    if arguments['--q9']:
        nbr_test = 100
        if arguments['<nbr_test>']:
            try:
                nbr_test = int(arguments['<nbr_test>'])
                if nbr_test <= 0:
                    raise ValueError
            except:
                print("nbr_test doit être un naturel non nul")
        try:
            a, alpha = float(arguments['<a>']), float(arguments['<alpha>'])
            print(q9(a, alpha, nbr_test))
        except:
            print("a et alpha et i doivent être des valeurs réelles")


    if arguments['--q10']:
        try:
            a = float(arguments['<a>'])
            i0 = float(arguments['<i0>'])
            find = float(arguments['<find>'])
            print(f"Valeur de alpha : {q10(a, i0, find)}")
        except:
            print("a et i et f doivent être des valeurs réelles")

    if arguments['--q11']:
        gamma = computeGamma(0.9, 0.1, 0.5, 0.3, (0.99, 0.01, 0), 0.15)
        T = computeTMax(0.9, 0.1, gamma, 0.5, 0.3, (0.99, 0.01, 0))
        print(f"Valeur de gamma : {gamma}")
        print(f"Valeur du maximum de T : {T}")
