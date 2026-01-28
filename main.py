from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
from plot import *
import numpy as np

base_table = RotTable("dna/table.json")
base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

nb_indiv = 100
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

res_exp_normal = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"roulette_exp_norm")
res_recuit_30 = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"roulette_exp", recuit=True)
best_exp_normal, score_exp_normal, _ = res_exp_normal
best_recuit30, score_recuit_30, _ = res_recuit_30
"""
score = fitness(best[-1].Rot_table,base_seq,nbcuts=0)
print(" score : ",best[-1].score,"score final : ",score," via type de selection : ","elitiste")
# traj_res = Traj3D(want_to_plot=True)
# traj_res.compute(base_seq,best[-1].Rot_table)
# traj_res.draw()
plot_with_slider(get_trajectories(best,base_seq))
save_trajectory_gif(get_trajectories(best,base_seq))
"""
plt.plot(score_recuit_30)
plt.plot(score_exp_normal)
plt.show()