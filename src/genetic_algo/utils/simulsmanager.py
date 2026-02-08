from .resultsmanager import save_simulation_data,load_simulation_data
from genetic_algo.core.algogenetique import *
from .plot import plot_with_slider,get_trajectories,save_trajectory_gif
import os

def simul_and_save_results(save_filename,dna_seq,params):
    """
    Run genetic algorithm simulation and save results to file.
    
    Args:
        save_filename: Output file path for simulation results (.pkl)
        dna_seq: DNA sequence string to optimize
        params: Dictionary of genetic algorithm parameters
                (nb_individus, nb_generations, taux_selec, selection_type, etc.)
    
    Automatically locates rotation table and saves best/worst scores per generation.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_src = os.path.dirname(os.path.dirname(current_dir)) # goes up to 'src'
    table_path = os.path.join(project_src, "genetic_algo", "dna", "table.json")
    res = AlgoGenetique(table_path,dna_seq,**params)
    print("Simulation Terminée ... Sauvegarde en cours")
    best,bscore,wscore = res
    save_simulation_data(save_filename,best,bscore,wscore,params,dna_seq)

def print_final_score_result(filename,dna_seq):
    """
    Load simulation results and print final fitness scores.
    
    Args:
        filename: Path to simulation pickle file
        dna_seq: DNA sequence string (for validation)
    
    Prints:
        - Best and worst fitness scores from last generation
        - Euclidean distance score for best individual
        - Simulation parameters used
    """
    res = load_simulation_data(filename,dna_seq)
    best,bscore,wscore,params =res
    eucl_score_fin = fitness(best[-1].Rot_table,dna_seq,nbcuts=0)
    print(f"Selon fitness, meilleur score final: {bscore[-1]}, pire score final: {wscore[-1]} \n meilleur score final en distance euclidienne : {eucl_score_fin}  Avec les paramètres : {params}")

def load_and_visualise_timeline(filename,dna_seq):
    """
    Load simulation and display interactive 3D trajectory evolution.
    
    Args:
        filename: Path to simulation pickle file
        dna_seq: DNA sequence string (for validation)
    
    Opens interactive plot with slider to navigate through generations.
    """
    res = load_simulation_data(filename,dna_seq)
    best,_,_,_ = res
    print("Chargement des données terminé, création du gif en cours..")
    plot_with_slider(get_trajectories(best,dna_seq))

def load_and_save_gif(save_filename,load_filename,dna_seq,fps=10):
    """
    Load simulation and export trajectory evolution as animated GIF.
    
    Args:
        save_filename: Output GIF file path
        load_filename: Path to simulation pickle file
        dna_seq: DNA sequence string (for validation)
        fps: Frames per second for animation (default: 10)
    """
    res = load_simulation_data(load_filename,dna_seq)
    best,_,_,_ = res
    print("Chargement des données terminé, création du gif en cours..")
    save_trajectory_gif(get_trajectories(best,dna_seq),save_filename,fps=fps)