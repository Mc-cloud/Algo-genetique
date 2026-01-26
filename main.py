from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique

base_table = RotTable("dna/table.json")
base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) #exemple utilisé de dinucléotide

nb_indiv = 10
nb_generations = 10
taux_selec = 0.5

res = AlgoGenetique("dna/table.json",base_seq,nb_indiv,nb_generations,taux_selec,"elitiste")
traj_res = Traj3D(want_to_plot=True)
traj_res.compute(base_seq,res.Rot_table)