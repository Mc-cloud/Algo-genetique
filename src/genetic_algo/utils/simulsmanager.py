from .resultsmanager import save_simulation_data,load_simulation_data
from genetic_algo.core.algogenetique import *
from .plot import plot_with_slider,get_trajectories,save_trajectory_gif
def simul_and_save_results(save_filename,dna_seq,params):
    """
    Permet de simuler selon les paramètres, et l'entrée (la séquence adn),
    l'algorithme génétique, et de stocker ses résultats dans le fichier correspondant (=save_filename)
    """
    res = AlgoGenetique("dna/table.json",dna_seq,**params)
    print("Simulation Terminée ... Sauvegarde en cours")
    best,bscore,wscore = res
    save_simulation_data(save_filename,best,bscore,wscore,params,dna_seq)

def print_final_score_result(filename,dna_seq):
    res = load_simulation_data(filename,dna_seq)
    best,bscore,wscore,params =res
    eucl_score_fin = fitness(best[-1].Rot_table,dna_seq,nbcuts=0)
    print(f"Selon fitness, meilleur score final: {bscore[-1]}, pire score final: {wscore[-1]} \n meilleur score final en distance euclidienne : {eucl_score_fin}  Avec les paramètres : {params}")

def load_and_visualise_timeline(filename,dna_seq):
    res = load_simulation_data(filename,dna_seq)
    best,_,_,_ = res
    print("Chargement des données terminé, création du gif en cours..")
    plot_with_slider(get_trajectories(best,dna_seq))

def load_and_save_gif(save_filename,load_filename,dna_seq,fps=10):
    res = load_simulation_data(load_filename,dna_seq)
    best,_,_,_ = res
    print("Chargement des données terminé, création du gif en cours..")
    save_trajectory_gif(get_trajectories(best,dna_seq),save_filename,fps=fps)