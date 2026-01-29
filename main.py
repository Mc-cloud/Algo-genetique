from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
from plot import *
import numpy as np
from simulsmanager import *
from benchmark import grid_search_params_save,grid_search_compare

base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide
name = "data_algo/gsearch1"


# Les différentes configurations de paramètres à simuler, ou comparer
params_to_search = {
                        "nb_individus":[150],
                        "nb_generations":[20],
                        "taux_selec":[0.3],
                        "selection_type":["elitiste","tournament","rang_geo","roulette_exp_norm"],
                        "poisson":[True,False],
                        "nb_cuts":[0],
                        "nb_append":[1],
                        "recuit":[False]
                                    }

#Simule les différentes configurations et sauvegarde les résultats
#grid_search_params_save(name,base_seq,params_to_search)

#Compare les différentes configurations déjà simulées et sauvegardées, et affiche les résultats de la meilleure configuration
grid_search_compare(name,base_seq,params_to_search)

