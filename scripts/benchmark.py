import os
import sys
import time
import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# Add project root to path to ensure imports work
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Imports
from genetic_algo.core.algogenetique import AlgoGenetique, generate_pop
from genetic_algo.utils.simulsmanager import simul_and_save_results
from genetic_algo.utils.resultsmanager import load_simulation_data
from genetic_algo.dna.Traj3D import Traj3D
from genetic_algo.utils.plot import get_indicators, plot_with_slider, get_trajectories, plot_three_indicators, save_trajectory_gif



def load_fasta(file_path):
    """Read FASTA file and return DNA sequence (headers excluded)."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return "".join([line.strip() for line in lines if not line.startswith(">")])

def run_comparison_benchmark(dna_seq, table_path):
    """
    Compare performance of different selection strategies.
    
    Tests: elitist, tournament, linear roulette, and exponential roulette selection.
    Returns list of tuples: (strategy_name, best_score, duration, best_individual_list)
    """
    print(f"\n{'='*60}")
    print(f"🚀 RUNNING COMPARISON BENCHMARK")
    print(f"DNA Sequence Length: {len(dna_seq)} bases")
    print(f"{'='*60}\n")

    strategies = ["elitiste", "tournament", "roulette_lin", "roulette_exp"]
    results = []

    params = {
        "filename": table_path,
        "dna_seq": dna_seq,
        "nb_individus": 1000,
        "nb_generations": 200,   
        "taux_selec": 0.15,
        "poisson": False,
        "nb_cuts": 0,
        "nb_append": 1,
    }

    master_pop = generate_pop(
        params["nb_individus"], 
        rot_table_path = table_path, 
        dna_seq = dna_seq, 
        nb_cuts = 0, 
        nb_append = 1
    )

    print(f"{'Strategy':<15} | {'Best Score':<30} | {'Time (s)':<10} | {'Status'}")
    print("-" * 55)

    for strat in strategies:
        start_time = time.time()
        
        best_list, best_scores, worst_scores = AlgoGenetique(
            filename = table_path,
            dna_seq = dna_seq,
            initial_population = master_pop,
            nb_individus = params['nb_individus'],
            taux_selec = params['taux_selec'],
            nb_cuts = params['nb_cuts'],
            nb_append = params['nb_append'],
            selection_type = strat,
            nb_generations = params['nb_generations']
        )
        
        duration = time.time() - start_time
        best_score = best_scores[-1]
        
        print(f"{strat:<15} | {best_score:<25.10e} | {duration:<10.2f} | ✅ Done")
        results.append((strat, best_score, duration, best_list))

    # Summary
    print("-" * 55)
    best_run = min(results, key=lambda x: x[1])
    print(f"\n🏆 WINNER: {best_run[0]} (Score: {best_run[1]:.10e})")

    return results

def run_grid_search(dna_seq, table_path, base_save_filename):
    """
    Exhaustive hyperparameter grid search.
    
    Tests all combinations of population sizes, selection rates, strategies,
    cuts, and appends. Evaluates three trajectory metrics:
    - Distance: Euclidean distance between start and end (closure quality)
    - Norm difference: ||v_start - v_end|| (continuity)
    - Dot product: v_start · v_end (alignment)
    
    Outputs: 'Evolution_metrique_genetique.png' with three metric plots.
    """
    print(f"\n{'='*60}")
    print(f"🧪 RUNNING GRID SEARCH (This may take a while...)")
    print(f"{'='*60}\n")

    # Define the search space
    params_listed = {
        "nb_individus": [100, 200, 500],
        "nb_generations": [150],
        "taux_selec": [0.15, 0.4, 0.5],
        "selection_type": ["elitiste", "roulette_exp", "rang_geo", "tournoi_elitiste"],
        "poisson": [False],
        "nb_cuts": [0, 1, 2],
        "nb_append": [1, 2],
        "recuit": [False]
    }

    master_pop = generate_pop(max(params_listed["nb_individus"]), 
        rot_table_path = table_path, 
        dna_seq = dna_seq, 
        nb_cuts = 0, 
        nb_append = 1
        )

    total_sims = np.prod([len(v) for v in params_listed.values()])
    print(f"Estimated number of simulations: {total_sims}")

    histories = {}

    count = 0

    traj_tool = Traj3D()

    eval_seq = dna_seq + dna_seq[0] + dna_seq[1]

    for nb_ind in params_listed["nb_individus"]:
        pop_init = master_pop[:nb_ind]
        for nb_gen in params_listed["nb_generations"]:
            for rate in params_listed["taux_selec"]:
                for sel_type in params_listed["selection_type"]:
                    for is_recuit in params_listed["recuit"]:
                        for cut in params_listed["nb_cuts"]:
                            for append in params_listed["nb_append"]:
                                config_key = (nb_ind , nb_gen, rate, sel_type, is_recuit, cut, append)
                                count += 1
                                print(f"[{count}/{total_sims}] Simulating: {sel_type}, Pop: {nb_ind}, Recuit: {is_recuit}, Number of cuts : {cut}, Number of appends : {append}")
                                
                                curr_params = {
                                    "nb_individus": nb_ind,
                                    "nb_generations": nb_gen,
                                    "taux_selec": rate,
                                    "selection_type": sel_type,
                                    "poisson": False,
                                    "nb_cuts": cut,
                                    "nb_append": append,
                                    "recuit": is_recuit
                                }

                                res = AlgoGenetique("src/genetic_algo/dna/table.json", dna_seq,initial_population = pop_init, **curr_params)
                                bests, best_scores, worst_scores = res

                                histories[config_key] = {'dist' : [], 'norm' : [], 'ps' : []}
                                
                                for best in bests : 
                                    traj_tool.compute(eval_seq, best.Rot_table)
                                    coords = traj_tool.getTraj()

                                    dist, norm_diff, dot_prod = get_indicators(coords)

                                    histories[config_key]['dist'].append(dist)
                                    histories[config_key]['norm'].append(norm_diff)
                                    histories[config_key]['ps'].append(dot_prod)

    sorted_histories = sorted(histories.items(), key = lambda item : item[1]['dist'][-1])

    top_10 = sorted_histories[:10]

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
                                

def main():
    """
    Genetic algorithm benchmark tool with three modes:
    
    compare: Quick strategy comparison (default)
    grid: Exhaustive hyperparameter search
    plot: Visualize existing results
    
    Usage:
        python benchmark.py --mode compare
        python benchmark.py --mode grid --file data/raw/my_sequence.fasta
    """

    parser = argparse.ArgumentParser(description="Genetic Algorithm Benchmark Tool")
    parser.add_argument('--mode', type=str, choices=['compare', 'grid', 'plot'], default='compare', 
                        help='Mode: "compare" (speed test), "grid" (exhaustive search), "plot" (visualize existing results)')
    parser.add_argument('--file', type=str, default='data/raw/plasmid_8k.fasta', 
                        help='Path to the fasta file')
    parser.add_argument('--table', type=str, default='src/genetic_algo/dna/table.json', 
                        help='Path to the rotation table JSON')
    
    args = parser.parse_args()

    # 1. Load Data
    full_path_fasta = os.path.join(project_root, args.file)
    full_path_table = os.path.join(project_root, args.table)
    
    if not os.path.exists(full_path_fasta):
        print(f"❌ Error: File not found at {full_path_fasta}")
        return

    dna_seq = load_fasta(full_path_fasta)

    # 2. Execute Mode
    if args.mode == 'compare':
        run_comparison_benchmark(dna_seq, full_path_table)
        
    elif args.mode == 'grid':
        output_base = os.path.join(project_root, 'data', 'processed', 'benchmark_')
        run_grid_search(dna_seq, full_path_table, output_base)
        
    elif args.mode == 'plot':
        print("To use plot mode, point to a specific saved simulation folder in the code.")
        # You can expand this part to load a specific file if needed.

if __name__ == "__main__":
    main()
