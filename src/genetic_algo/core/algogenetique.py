import src.genetic_algo.dna.RotTable as RotTable
import src.genetic_algo.dna.Traj3D as Traj3D
import numpy as np
import random
from src.genetic_algo.core.fitness import fitness,fitness_basic
from src.genetic_algo.core.selection import selection
import copy
from json import load as json_load

str_data = None
Rot_data = None
nb_cut = None
nbappend = None
beta = 1
big_mut = 2
class Individu:
    def __init__(self, Table_rot):
        """
        Initialise un individu avec une table de rotations.
        Args:
            Table_rot: Dictionnaire contenant les paramètres de rotation pour chaque dinucléotide
        """
        self.Rot_table = RotTable.RotTable(Table_rot)
        self.score = self.fit()


    def __add__(self, other): #fonction permettant d'acoupler deux individus
        """
        Opérateur de croisement (crossover) : combine deux individus parents.
        Effectue un croisement pondéré par les scores des parents. Le parent avec le
        meilleur score (plus bas) a plus d'influence sur la descendance.
        Args:
            other: Autre individu parent
        Returns:
            Nouvel individu enfant issu du croisement
        Note:
            - alpha représente la contribution du parent 'other'
            - beta contrôle l'intensité du mélange génétique
            - Chaque paramètre a 50% de chance de provenir majoritairement de chaque parent
        """
        s_score,o_score = self.score,other.score
        alpha = 0.5
        if s_score >0 or o_score >0 :
            alpha = (o_score)/(s_score+o_score) 
        Table ={}
        for XY in self.Rot_table.rot_table:
            Ls = self.Rot_table.rot_table[XY].copy()
            Lo = other.Rot_table.rot_table[XY].copy()
            for i in range(3):
                if random.random() < alpha :
                    Ls[i] = (alpha*Ls[i] +(1-alpha)*Lo[i])*beta + (1-beta)* Ls[i] 
                else :
                    Ls[i] = (alpha*Ls[i] +(1-alpha)*Lo[i])*beta + (1-beta)* Lo[i] 
            Table[XY] = Ls
        return Individu(Table)
    
    def mutation(self,mutrate,sigma):
        #global Rot_data
        if 0<=mutrate <=1:
            Table =copy.deepcopy(self.Rot_table.rot_table)
            for XY in Table:
                L = Table[XY]
                for i in range(3):
                    if random.random() < mutrate :  #petites mutations fréquentes
                        tamp = L[i] + random.normalvariate(0,sigma*L[i+3])
                        L[i] = max(Rot_data[XY][i]-Rot_data[XY][i+3],min(tamp,Rot_data[XY][i]+Rot_data[XY][i+3]))
                    if random.random() < mutrate/big_mut : # grosses mutations rares
                        tamp = L[i] + random.normalvariate(0,sigma*np.sqrt(np.sqrt(big_mut))*L[i+3])
                        L[i] = max(Rot_data[XY][i]-Rot_data[XY][i+3],min(tamp,Rot_data[XY][i]+Rot_data[XY][i+3]))
                Table[XY] = L
            self.Rot_table = RotTable.RotTable(Table)
            self.score = self.fit()

    def fit(self) -> float:
        """
        Calcule le score de fitness de l'individu.
        """
        return fitness(self.Rot_table,str_data,nbcuts=nb_cut,nbappend=nbappend)

    def __lt__(self,other):
        """
        Opérateur de comparaison : permet de trier les individus par score.
        """
        return self.score<other.score
    

def AlgoGenetique(filename : str,dna_seq: str, nb_individus,nb_generations,taux_selec,selection_type : str,poisson=False,nb_cuts = 0,nb_append = 1,recuit=False,beta_reproduction = 0.7,mutrate = 0.02, big_mutation = 20) :
    """
    Algorithme génétique pour optimiser les paramètres de rotation de l'ADN.
    
    Cherche à trouver une table de rotations qui minimise le score de fermeture
    de la structure ADN (pour former des plasmides circulaires stables).
    
    Args:
        filename: Chemin vers le fichier JSON contenant la table de rotations initiale
        dna_seq: Séquence ADN à optimiser
        nb_individus: Taille de la population
        nb_generations: Nombre de générations à simuler
        taux_selec: Taux de sélection (proportion d'individus gardés comme géniteurs)
        selection_type: Type de sélection ("elitiste", "tournament", "roulette", etc.)
        poisson: Si True, le nombre de géniteurs suit une loi de Poisson (défaut: False)
        nb_cuts: Nombre de points de coupure pour le calcul de fitness (défaut: 0)
        nb_append: Nombre de nœuds à ajouter pour tester la fermeture (défaut: 1)
        recuit: Si True, active le recuit simulé (température décroissante) (défaut: False)
        beta_reproduction: Coefficient de mélange lors du croisement (défaut: 0.7)
        mutrate: Taux de mutation initial (défaut: 0.02)
        big_mutation: Facteur pour les grosses mutations (défaut: 20)
    
    Returns:
        Tuple contenant :
        - Liste des meilleurs individus à chaque génération
        - Liste des scores du meilleur individu à chaque génération
        - Liste des scores du pire individu à chaque génération
    """
    rot_table = json_load(open(filename))
    global str_data,Rot_data, nb_cut, nbappend, big_mut, beta
    str_data = dna_seq
    Rot_data = rot_table
    nb_cut = nb_cuts
    nbappend = nb_append 
    beta = beta_reproduction
    big_mut = big_mutation
    def New_pop():
        def New_individu():
            Table_rot =  copy.deepcopy(rot_table)
            for XY in Table_rot:
                L = Table_rot[XY]
                for i in range(3):
                        L[i] += random.uniform(-Rot_data[XY][i+3],Rot_data[XY][i+3])
                Table_rot[XY] = L
            return Individu(Table_rot)
        L = [New_individu() for _ in range (nb_individus)]
        return L
    
    Population = New_pop()
    Best_indiv_list = [min(Population)]
    Best_indiv_score_list = [Best_indiv_list[0].score]
    worst_indiv_score_list = [max(Population).score]

    for i in range(nb_generations):
        print("itération :", i+1, "/", nb_generations)
        if poisson: # si l'on s'est mit en mode processur aléatoire de poisson, on aura une nombre de géniteurs suivant une loi de Poisson(taux_selec*nb_indiv) (on prendra le max avec 2, pour avoir assez de géniteurs)
            if recuit:
                Geniteurs = selection(Population,max(2,np.random.poisson(taux_selec*nb_individus))/nb_individus,selection_type, n=i)
            else:
                Geniteurs = selection(Population,max(2,np.random.poisson(taux_selec*nb_individus))/nb_individus,selection_type)
        else:
            if recuit:
                Geniteurs = selection(Population,taux_selec,selection_type, n=i)
            else:
                Geniteurs = selection(Population,taux_selec,selection_type)
        A =[Geniteurs[k].score for k in range(len(Geniteurs))]
        print("fit : ", np.sum(A))
        Population = Geniteurs.copy()
        while len(Population) < nb_individus:
            individu = random.choice(Geniteurs)+random.choice(Geniteurs)
            individu.mutation(mutrate*(1-i/nb_generations),(1-i/nb_generations)*0.5)  
            Population.append(individu)
        best_indiv = min(Population)
        worst_indiv = max(Population)
        print(f"Meilleur pour iter {i+1} : {best_indiv.score}")
        print(f"Pire pour iter {i+1} : {worst_indiv.score}")
        Best_indiv_list.append(best_indiv)
        Best_indiv_score_list.append(best_indiv.score)
        worst_indiv_score_list.append(worst_indiv.score)

    return Best_indiv_list,Best_indiv_score_list,worst_indiv_score_list