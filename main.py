from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
import numpy as np

base_table = RotTable("dna/table.json")
base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

nb_indiv = 300
nb_generations = 25
taux_selec = 0.5


'''
for selection_type in selections_dic.keys():
    results = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,selection_type,poisson=False)
    res = results[-1]
    score = fitness(res.Rot_table,base_seq,nbcuts=0)
    print(" score : ",res.score,"score final : ",score," via type de selection : ",selection_type)
#     # traj_res = Traj3D(want_to_plot=True)
#     # traj_res.compute(base_seq,res.Rot_table)
#     # traj_res.draw() #'''

res,_,_ = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"elitiste")
score = fitness(res[-1].Rot_table,base_seq,nbcuts=0)
print(" score : ",res[-1].score,"score final : ",score," via type de selection : ","elitiste")
traj_res = Traj3D(want_to_plot=True)
traj_res.compute(base_seq,res[-1].Rot_table)
traj_res.draw()
