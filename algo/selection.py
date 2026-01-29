import random
import numpy as np

def selection_elitiste(list_ind, taux_selec):
    """
    Sélection élitiste : conserve les meilleurs individus.
    
    Args:
        list_ind: Liste des individus à sélectionner
        taux_selec: Taux de sélection (proportion d'individus à garder)
    
    Returns:
        Liste des meilleurs individus selon le taux de sélection
    """
    list_ind_sort = list(list_ind)
    list_ind_sort.sort()

    n = len(list_ind)

    selection = list_ind_sort[:int(n*taux_selec)]

    return selection

def selection_tournament_elitiste(list_ind, taux_selec, p = 0.001):
    """
    Sélection par tournoi : compare deux individus aléatoires et sélectionne le meilleur.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
        p: Probabilité de sélectionner le moins bon individu (diversité)
    
    Returns:
        Liste d'individus sélectionnés par tournoi + 10% des meilleurs
    """
    selection = []
    n_list = list(list_ind)
    n = len(n_list)
    k = int(n*taux_selec)

    for _ in range(k):
        q = random.random()
        x = random.sample(n_list, k=2)
        if (x[0].score > x[1].score and q > p) or (x[0].score < x[1].score and q <= p) :
            selection.append(x[1])
        elif  (x[0].score > x[1].score and q <= p) or (x[0].score < x[1].score and q > p): 
            selection.append(x[0])

    list_ind_sort = list(list_ind)
    list_ind_sort.sort()

    for ind in list_ind_sort[:int(n*0.1)]:
        if ind not in selection:
            selection.append(ind)

    return selection

def selection_tournament(list_ind, taux_selec, p = 0.001):
    """
    Sélection par tournoi : compare deux individus aléatoires et sélectionne le meilleur.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
        p: Probabilité de sélectionner le moins bon individu (diversité)
    
    Returns:
        Liste d'individus sélectionnés par tournoi + 10% des meilleurs
    """
    selection = []
    n_list = list(list_ind)
    n = len(n_list)
    k = int(n*taux_selec)

    for _ in range(k):
        q = random.random()
        x = random.sample(n_list, k=2)
        if (x[0].score > x[1].score and q > p) or (x[0].score < x[1].score and q <= p) :
            selection.append(x[1])
        elif  (x[0].score > x[1].score and q <= p) or (x[0].score < x[1].score and q > p): 
            selection.append(x[0])

    return selection

def selection_roulette(list_ind, taux_selec):
    """
    Sélection par roulette : probabilité inversement proportionnelle au score.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
    
    Returns:
        Liste d'individus sélectionnés selon leurs probabilités
    """
    total = sum([ind.score for ind in list_ind])
    proba = [1-ind.score/total for ind in list_ind]
    
    select = random.choices(list_ind, proba, k = int(len(list_ind)*taux_selec))

    return select

def selection_roulette_exp(list_ind, taux_selec, temp=1000):
    """
    Sélection par roulette avec fonction exponentielle (contrôle de la pression de sélection).
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
        temp: Température (contrôle la pression de sélection, plus elle est élevée, moins la sélection est sévère)
    
    Returns:
        Liste d'individus sélectionnés avec distribution exponentielle
    """
    min_score = min([ind.score for ind in list_ind])
    proba = [np.exp((min_score**2-ind.score**2)/temp) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)))
    return select

def selection_roulette_exp_normal(list_ind, taux_selec):
    """
    Sélection par roulette exponentielle normalisée (température = min_score).
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
    
    Returns:
        Liste d'individus sélectionnés avec distribution exponentielle normalisée
    """
    min_score = min([ind.score for ind in list_ind])
    proba = [np.exp((min_score**2-ind.score**2)/min_score) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)))
    return select

def selection_rang_reel(list_ind, taux_selec):
    """
    Sélection par rang linéaire : probabilité proportionnelle au rang.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
    
    Returns:
        Liste d'individus sélectionnés selon leur rang
    """
    list_ind.sort(reverse = True)
    
    proba = [i + 1 for i in range(len(list_ind))]
    
    n_parents = len(list_ind)
    select = random.choices(list_ind, weights=proba, k=int(taux_selec*n_parents))
    
    return select

def selection_rang_geometrique(list_ind, taux_selec, q=0.3):
    """
    Sélection par rang géométrique : probabilité décroissante géométriquement.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
        q: Paramètre de pression de sélection (0 < q < 1)
           Plus q est élevé, plus la sélection favorise les meilleurs
    
    Returns:
        Liste d'individus sélectionnés selon une distribution géométrique
    """
    list_ind.sort()

    n = len(list_ind)

    proba = [q * (1 - q)**i for i in range(n)]
    
    # 3. Select parents
    k = int(len(list_ind) * taux_selec)
    select = random.choices(list_ind, weights=proba, k=k)
    
    return select

selections_dic = {"tournament_elitiste" : selection_tournament_elitiste, "elitiste":selection_elitiste,"tournament":selection_tournament,"roulette":selection_roulette,"rang_reel":selection_rang_reel,"rang_geo":selection_rang_geometrique, "roulette_exp" : selection_roulette_exp, "roulette_exp_norm": selection_roulette_exp_normal}
def selection(list_ind,taux_selec,select_type, n=None):
    """
    Fonction générique de sélection : appelle la méthode appropriée.
    
    Args:
        list_ind: Liste des individus
        taux_selec: Taux de sélection
        select_type: Type de sélection (clé du dictionnaire selections_dic)
        n: Numéro de génération (optionnel, utilisé pour ajuster la température)
    
    Returns:
        Liste d'individus sélectionnés selon la méthode choisie
    """
    if n:
        if select_type == "roulette_exp":
            return selections_dic[select_type](list_ind,taux_selec, temp = max(100,700-n*30))
    if select_type in selections_dic:
        return selections_dic[select_type](list_ind,taux_selec)
    return selection_elitiste(list_ind,taux_selec)



