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

params = {
    "nb_individus":8,
    "nb_generations":10,
    "taux_selec":0.5,
    "selection_type":"elitiste",
    "poisson":False,
    "nb_cuts":0,
    "nb_append":1,
    "recuit":False
    }



#simul_and_save_results("data_algo/exemple1",base_seq,params)

#load_and_visualise_timeline("data_algo/benchmark_nbcuts0_nbappend1",base_seq)

#load_and_save_gif("gifs/benchmark_nbcuts0_nbappend1.gif","data_algo/benchmark_nbcuts0_nbappend1",base_seq)

name = "data_algo/gsearch2"
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
grid_search_params_save(name,base_seq,params_to_search)
grid_search_compare(name,base_seq,params_to_search)

