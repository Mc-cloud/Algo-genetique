import random
import numpy as np

def selection_elitiste(list_ind, taux_selec):
    """
    Elitist selection: keeps the best individuals.
    
    Args:
        list_ind: List of individuals to select
        selection_rate: Selection rate (proportion of individuals to keep)
    
    Returns:
        List of the best individuals according to the selection rate
    """
    list_ind_sort = list(list_ind)
    list_ind_sort.sort()

    n = len(list_ind)

    selection = list_ind_sort[:int(n*taux_selec)]

    return selection

def selection_tournament_elitiste(list_ind, taux_selec, p = 0.001):
    """
    Elitist Tournament selection: compares two random individuals and selects the best one.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
        p: Probability of selecting the worst individual (diversity)
    
    Returns:
        List of individuals selected by tournament + 10% of the best ones
    """
    selection = []
    n_list = list(list_ind)
    n = len(n_list)
    k = int(n*taux_selec)

    for _ in range(k):
        q = random.random()
        x = random.sample(n_list, k=2)
        if (x[0].score >= x[1].score and q > p) or (x[0].score < x[1].score and q <= p) :
            selection.append(x[1])
        elif  (x[0].score >= x[1].score and q <= p) or (x[0].score < x[1].score and q > p): 
            selection.append(x[0])

    list_ind_sort = list(list_ind)
    list_ind_sort.sort()

    for ind in list_ind_sort[:int(n*0.1)]:
        if ind not in selection:
            selection.append(ind)

    return selection

def selection_tournament(list_ind, taux_selec, p = 0.001):
    """
    Tournament selection: compares two random individuals and selects the best one.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
        p: Probability of selecting the worst individual (diversity)
    
    Returns:
        List of individuals selected by tournament
    """
    selection = []
    n_list = list(list_ind)
    n = len(n_list)
    k = int(n*taux_selec)

    for _ in range(k):
        q = random.random()
        x = random.sample(n_list, k=2)
        if (x[0].score >= x[1].score and q > p) or (x[0].score <= x[1].score and q <= p) :
            selection.append(x[1])
        elif  (x[0].score >= x[1].score and q <= p) or (x[0].score <= x[1].score and q > p): 
            selection.append(x[0])

    return selection

def selection_roulette(list_ind, taux_selec):
    """
    Selection by roulette wheel: probability inversely proportional to score.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
    
    Returns:
        List of individuals selected according to their probabilities
    """
    total = sum([ind.score for ind in list_ind])
    proba = [1-ind.score/total for ind in list_ind]
    
    select = random.choices(list_ind, proba, k = int(len(list_ind)*taux_selec))

    return select

def selection_roulette_exp(list_ind, taux_selec, temp=1000):
    """
    Selection by roulette wheel with exponential function (selection pressure control).
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
        temp: Temperature (controls selection pressure; the higher it is, the less severe the selection)
    
    Returns:
        List of selected individuals with exponential distribution
    """
    min_score = min([ind.score for ind in list_ind])
    proba = [np.exp((min_score**2-ind.score**2)/temp) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)))
    return select

def selection_roulette_exp_normal(list_ind, taux_selec):
    """
    Selection by normalized exponential roulette wheel (temperature = min_score).
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
    
    Returns:
        List of selected individuals with normalized exponential distribution
    """
    min_score = min([ind.score for ind in list_ind])
    proba = [np.exp((min_score**2-ind.score**2)/(2*min_score)) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)))
    return select

def selection_rang_reel(list_ind, taux_selec):
    """
    Linear rank selection: probability proportional to rank.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
    
    Returns:
        List of individuals selected according to their rank
    """
    list_ind.sort(reverse = True)
    
    proba = [i + 1 for i in range(len(list_ind))]
    
    n_parents = len(list_ind)
    select = random.choices(list_ind, weights=proba, k=int(taux_selec*n_parents))
    
    return select

def selection_rang_geometrique(list_ind, taux_selec, q=0.3):
    """
    Selection by geometric rank: geometrically decreasing probability.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
        q: Selection pressure parameter (0 < q < 1)
           The higher q is, the more selection favors the best individuals.
    
    Returns:
        List of individuals selected according to a geometric distribution
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
    Generic selection function: calls the appropriate method.
    
    Args:
        list_ind: List of individuals
        selection_rate: Selection rate
        select_type: Selection type (key from the selections_dic dictionary)
        n: Generation number (optional, used to adjust the temperature)
    
    Returns:
        List of individuals selected according to the chosen method
    """
    if n:
        if select_type == "roulette_exp":
            return selections_dic[select_type](list_ind,taux_selec, temp = max(100,700-n*30))
    if select_type in selections_dic:
        return selections_dic[select_type](list_ind,taux_selec)
    return selection_elitiste(list_ind,taux_selec)



