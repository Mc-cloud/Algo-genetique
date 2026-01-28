from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
from plot import *
import numpy as np
from simulsmanager import *

base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 63ff9f3f8ff91926f3d7d475de3b1a8b16d99330
params = {
    "nb_individus":100,
    "nb_generations":20,
    "taux_selec":0.5,
    "selection_type":"elitiste",
    "poisson":False,
<<<<<<< HEAD
    "nb_cuts":0,
    "nb_append":1
    }
=======
nb_indiv = 100
nb_generations = 20
taux_selec = 0.5
>>>>>>> 381d799954e0adaa6bccad1350495eb29f4bb1c2
=======
    "nb_cuts": 0,
    "nb_append" : 1
    }
>>>>>>> 63ff9f3f8ff91926f3d7d475de3b1a8b16d99330

#simul_and_save_results("data_algo/exemple1",base_seq,params)

#load_and_visualise_timeline("data_algo/exemple1",base_seq)

load_and_save_gif("gifs/exemple1","data_algo/exemple1",base_seq)
