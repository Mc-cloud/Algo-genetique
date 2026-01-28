from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
import numpy as np

base_table = RotTable("dna/table.json")
base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

nb_indiv = 150
nb_generations = 20
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
T = [(1,5),(1,1)]
L = []
for a,b in T :
    res = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"elitiste",nb_cuts = a,nb_append = b)
    bests, _, _ = res
    best = bests[-1]
    traj = Traj3D()
    traj.compute(base_seq+base_seq[0]+base_seq[1], best.Rot_table)
    coords = traj.getTraj()
    dist = np.linalg.norm(coords[0]-coords[-2])
    v1 = coords[-1]-coords[-2]
    v2 = coords[1] -coords[0]
    v1,v2 = v1/np.linalg.norm(v1),v2/np.linalg.norm(v2)
    print("dist :", dist, "norm:", np.linalg.norm(v1-v2), "ps :", np.dot(v2,v1))

'''
score = fitness(best.Rot_table,base_seq,nbcuts=0)
print(" score : ",best.score,"score final : ",score," via type de selection : ","elitiste")
traj_res = Traj3D(want_to_plot=True)
traj_res.compute(base_seq,best[-1].Rot_table)
traj_res.draw() #'''
