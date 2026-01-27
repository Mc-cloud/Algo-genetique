import random
import numpy as np

def selection_elitiste(list_ind, taux_selec):
    list_ind_sort = list(list_ind)
    list_ind_sort.sort()

    n = len(list_ind)

    selection = list_ind_sort[:int(n*taux_selec)]

    return selection

def selection_tournament(list_ind, taux_selec):
    selection = []
    n_list = list(list_ind)
    n = len(n_list)
    k = int(n*taux_selec)

    for _ in range(k):
        x = random.sample(n_list, k=2)
        if x[0].score > x[1].score:
            selection.append(x[1])
        else : 
            selection.append(x[0])
    return selection

def selection_roulette(list_ind, taux_selec):
    total = sum([ind.score for ind in list_ind])
    proba = [1-ind.score/total for ind in list_ind]
    
    select = random.choices(list_ind, proba, k = int(len(list_ind)*taux_selec))

    return select

def selection_roulette_exp(list_ind, taux_selec, temp=1000):
    min_score = min([ind.score for ind in list_ind])
    proba = [np.exp((min_score**2-ind.score**2)/temp) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)))
    return select

def selection_rang_reel(list_ind, taux_selec):
    list_ind.sort(reverse = True)
    
    proba = [i + 1 for i in range(len(list_ind))]
    
    n_parents = len(list_ind)
    select = random.choices(list_ind, weights=proba, k=int(taux_selec*n_parents))
    
    return select

def selection_rang_geometrique(list_ind, taux_selec, q=0.3):
    """
    q: pression de selec (0 < q < 1).
    """
    list_ind.sort()

    n = len(list_ind)

    proba = [q * (1 - q)**i for i in range(n)]
    
    # 3. Select parents
    k = int(len(list_ind) * taux_selec)
    select = random.choices(list_ind, weights=proba, k=k)
    
    return select

selections_dic = {"elitiste":selection_elitiste,"tournament":selection_tournament,"roulette":selection_roulette,"rang_reel":selection_rang_reel,"rang_geo":selection_rang_geometrique, "roulette_exp" : selection_roulette_exp}
def selection(list_ind,taux_selec,select_type, n=None):
    if n:
        if select_type == "roulette_exp":
            return selections_dic[select_type](list_ind,taux_selec, temp = max(100,1000-n*10))
    if select_type in selections_dic:
        return selections_dic[select_type](list_ind,taux_selec)
    return selection_elitiste(list_ind,taux_selec)



