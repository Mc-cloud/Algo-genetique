import genetic_algo.dna.RotTable as RotTable
import genetic_algo.dna.Traj3D as Traj3D
import numpy as np
import random
from genetic_algo.core.fitness import fitness,fitness_basic
from genetic_algo.core.selection import selection
import copy
from json import load as json_load


# =============================================================================
# GLOBAL VARIABLES
# =============================================================================

str_data = None
Rot_data = None
nb_cut = None
nbappend = None
beta = 1
big_mut = 2



class Individu:
    """
    Individual solution candidate in the genetic algorithm.
    Contains rotation parameters for all 16 dinucleotides (AA, AC, ..., TT).
    """
    def __init__(self, Table_rot): 
        """Initialize individual with rotation table and compute fitness score."""
        self.Rot_table = RotTable.RotTable(Table_rot)
        self.score = self.fit()


    def __add__(self, other): #fonction permettant d'acoupler deux individus
        """
        Crossover operator: combine two parents to create offspring.
        Uses weighted mixing based on parent fitness scores.
        Better parents (lower score) contribute more to offspring.
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
        """
        Apply mutations to rotation parameters.
        Two types: frequent small mutations and rare large mutations.
        Parameters are bounded within valid ranges from reference table.
        """

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
        """Calculate fitness score (lower is better)."""
        return fitness(self.Rot_table,str_data,nbcuts=nb_cut,nbappend=nbappend)

    def __lt__(self,other):
        """Enable sorting by fitness score."""
        return self.score<other.score
    
def generate_pop(nb_individus, rot_table_path, dna_seq, nb_cuts = 0, nb_append = 1):
    """
    Generate initial population with randomized rotation parameters.
    Each parameter is randomly perturbed within ±std_dev of reference value.
    """

    rot_table = json_load(open(rot_table_path))
    global str_data, Rot_data, nb_cut, nbappend
    str_data = dna_seq
    Rot_data = rot_table
    nb_cut = nb_cuts
    nbappend = nb_append

    def New_individu():
        Table_rot = copy.deepcopy(rot_table)

        for XY in Table_rot:
            L = Table_rot[XY]
            for i in range(3) : 
                L[i] += random.uniform(-Rot_data[XY][i+3], Rot_data[XY][i+3])
            Table_rot[XY] = L
        return Individu(Table_rot)
    return [New_individu() for _ in range(nb_individus)]
    

def AlgoGenetique(filename : str,dna_seq: str, 
                nb_individus,nb_generations,taux_selec,selection_type : str,
                poisson=False,nb_cuts = 0,nb_append = 1,recuit=False,beta_reproduction = 0.7,
                mutrate = 0.02, big_mutation = 20, initial_population = None) :
    """
    Genetic algorithm for optimizing DNA rotation parameters.
    
    Evolves population to minimize closure error of circular DNA structures.
    Process: selection → crossover → mutation → replacement.
    
    Args:
        filename: Path to reference rotation table JSON
        dna_seq: DNA sequence to optimize
        nb_individus: Population size
        nb_generations: Number of evolution iterations
        taux_selec: Selection rate (fraction kept as parents, 0-1)
        selection_type: "elitiste", "tournament", "roulette_lin", "roulette_exp", etc.
        poisson: Use Poisson distribution for parent count (default: False)
        nb_cuts: Cut points for trajectory calculation (default: 0)
        nb_append: Bases to append for closure test (default: 1)
        recuit: Enable simulated annealing - adaptive mutation (default: False)
        beta_reproduction: Crossover mixing coefficient (default: 0.7)
        mutrate: Initial mutation rate (default: 0.02)
        big_mutation: Large mutation factor (default: 20)
        initial_population: Pre-generated population (default: None, generates random)
    
    Returns:
        tuple: (best_individuals_per_gen, best_scores, worst_scores)
    """
    rot_table = json_load(open(filename))
    global str_data,Rot_data, nb_cut, nbappend, big_mut, beta
    str_data = dna_seq
    Rot_data = rot_table
    nb_cut = nb_cuts
    nbappend = nb_append 
    beta = beta_reproduction
    big_mut = big_mutation

    if initial_population is not None : 
        Population = copy.deepcopy(initial_population)
        for ind in Population:
            ind.score = ind.fit()
    else:
        Population = generate_pop(nb_individus, filename, dna_seq, nb_cuts, nb_append)

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