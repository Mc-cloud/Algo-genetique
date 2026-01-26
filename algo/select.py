import random
import numpy as np

def selection_eletiste(list_ind, taux_selec):
    list_ind_sort = sorted(list_ind)

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
    proba = []
    for ind in list_ind:
        prob = 1 - ind.score/total
        proba.append(prob)
    
    select = random.choices(list_ind, proba, k = int(len(list_ind)*taux_selec))

    return select

def selection_rang(list_ind, taux_selec):
    proba = [np.exp(-ind.score) for ind in list_ind]
    select = random.choices(list_ind, proba, k = int(taux_selec*len(list_ind)//2))
    return select

def selection_rang_reel(list_ind, taux_selec):
    list_ind_sort = sorted(list_ind, reverse=True)
    
    proba = [i + 1 for i in range(len(list_ind_sort))]
    
    n_parents = len(list_ind) // 2
    select = random.choices(list_ind_sort, weights=proba, k=int(taux_selec*n_parents))
    
    return select






