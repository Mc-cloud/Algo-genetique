import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)


from genetic_algo.dna.RotTable import *
from genetic_algo.dna.Traj3D import *
from genetic_algo.core.algogenetique import AlgoGenetique
from genetic_algo.core.selection import selections_dic
from genetic_algo.core.fitness import fitness
import numpy as np
import matplotlib.pyplot as plt

base_table = RotTable("src/genetic_algo/dna/table.json")
base_seq = ''.join([line.rstrip('\n') for line in open("data/raw/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

nb_indiv = 1000
nb_generations = 30
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
T = [(1,1), (0,1), (2,2), (2,1), (1,2)]
L = []

histories = {}

for a,b in T :
    res = AlgoGenetique("src/genetic_algo/dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"elitiste",nb_cuts = a,nb_append = b)
    bests, _, _ = res

    config_key = f"cuts = {a}, append = {b}"
    histories[config_key] = {'dist' : [], 'norm' : [], 'ps' : []}

    for best in bests : 
        traj = Traj3D()
        traj.compute(base_seq+base_seq[0]+base_seq[1], best.Rot_table)
        coords = traj.getTraj()
        dist = np.linalg.norm(coords[0]-coords[-2])
        v_end = coords[-1]-coords[-2]
        v_start = coords[1] -coords[0]
        v_end,v_start = v_end/np.linalg.norm(v_end),v_start/np.linalg.norm(v_start)
        norm_diff = np.linalg.norm(v_start - v_end)
        dot_prod = np.dot(v_start, v_end)

        histories[config_key]['dist'].append(dist)
        histories[config_key]['norm'].append(norm_diff)
        histories[config_key]['ps'].append(dot_prod)


fig, axs = plt.subplots(3,1, figsize = (10,12), sharex = True)

for label, data in histories.items():
    generations = range(len(data['dist']))

    axs[0].plot(generations, data['dist'], label=label)
    axs[0].set_ylabel('Distance')
    axs[0].set_title('Évolution de la distance de fermeture')
    axs[0].grid(True, linestyle='--')
    axs[0].legend()

    axs[1].plot(generations, data['norm'], label=label)
    axs[1].set_ylabel('$|\\vec{v}_{start} - \\vec{v}_{end}|$')
    axs[1].set_title('Continuité : Norme de la différence des directions')
    axs[1].grid(True, linestyle='--')

    axs[2].plot(generations, data['ps'], label=label)
    axs[2].set_ylabel('$\\vec{v}_{start} \\cdot \\vec{v}_{end}$')
    axs[2].set_xlabel('Génération')
    axs[2].set_title('Alignement : Produit scalaire des directions')
    axs[2].grid(True, linestyle='--')

plt.tight_layout()
plt.savefig('Evolution_metrique_genetique.png')
