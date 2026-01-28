"""
Compare fitness evolution across different selection types
"""

from dna.RotTable import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
import matplotlib.pyplot as plt
import numpy as np

import os

# Create the directory relative to where your script is
output_dir = "data_algo"
os.makedirs(output_dir, exist_ok=True)

# Load DNA sequence
base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:])

# Fixed parameters
nb_indiv = 150
nb_generations = 20
taux_selec = 0.5
nb_cuts = 0
nb_append = 1

# Get all available selection types
selection_types = list(selections_dic.keys())
print(f"Available selection types: {selection_types}\n")

# Store results
results_by_selection = {}

# Run genetic algorithm for each selection type
for selection_type in selection_types:
    print(f"\n{'='*70}")
    print(f"Running with selection type: {selection_type}")
    print(f"{'='*70}")
    
    bests, best_scores, worst_scores = AlgoGenetique(
        "dna/table.json",
        base_seq,
        nb_indiv,
        nb_generations,
        taux_selec,
        selection_type,
        poisson=False,
        nb_cuts=nb_cuts,
        nb_append=nb_append
    )
    
    results_by_selection[selection_type] = {
        'best_scores': best_scores,
        'worst_scores': worst_scores,
        'best_individuals': bests
    }
    
    print(f"\nFinal best score: {best_scores[-1]:.6f}")
    print(f"Final worst score: {worst_scores[-1]:.6f}")

# ============================================================================
# PLOT 1: Best fitness evolution comparison
# ============================================================================
plt.figure(figsize=(14, 8))
colors = plt.cm.tab10(np.linspace(0, 1, len(selection_types)))

for idx, selection_type in enumerate(selection_types):
    best_scores = results_by_selection[selection_type]['best_scores']
    generations = np.arange(len(best_scores))
    plt.plot(generations, best_scores, marker='o', linewidth=2.5, 
             markersize=5, label=selection_type, color=colors[idx], alpha=0.85)

plt.xlabel('Generation', fontsize=14, fontweight='bold')
plt.ylabel('Fitness Score (Best Individual)', fontsize=14, fontweight='bold')
plt.title('Fitness Evolution: Comparison of Selection Types', fontsize=16, fontweight='bold', pad=20)
plt.legend(fontsize=11, loc='best', framealpha=0.95, edgecolor='black')
plt.grid(True, alpha=0.3, linestyle='--', linewidth=1)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selection_comparison_best.png'), dpi=300, bbox_inches='tight')
print("\n✓ Plot saved: selection_comparison_best.png")

# ============================================================================
# PLOT 2: Best AND Worst fitness evolution (subplots)
# ============================================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 11))

# Top plot: Best fitness
for idx, selection_type in enumerate(selection_types):
    best_scores = results_by_selection[selection_type]['best_scores']
    generations = np.arange(len(best_scores))
    ax1.plot(generations, best_scores, marker='o', linewidth=2.5, 
             markersize=5, label=selection_type, color=colors[idx], alpha=0.85)

ax1.set_xlabel('Generation', fontsize=13, fontweight='bold')
ax1.set_ylabel('Fitness Score (Best)', fontsize=13, fontweight='bold')
ax1.set_title('Best Individual Fitness Evolution', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10, loc='best', framealpha=0.95, edgecolor='black')
ax1.grid(True, alpha=0.3, linestyle='--')

# Bottom plot: Worst fitness
for idx, selection_type in enumerate(selection_types):
    worst_scores = results_by_selection[selection_type]['worst_scores']
    generations = np.arange(len(worst_scores))
    ax2.plot(generations, worst_scores, marker='s', linewidth=2, 
             markersize=4, label=selection_type, color=colors[idx], alpha=0.85)

ax2.set_xlabel('Generation', fontsize=13, fontweight='bold')
ax2.set_ylabel('Fitness Score (Worst)', fontsize=13, fontweight='bold')
ax2.set_title('Worst Individual Fitness Evolution', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10, loc='best', framealpha=0.95, edgecolor='black')
ax2.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selection_comparison_combined.png'), dpi=300, bbox_inches='tight')
print("✓ Plot saved: selection_comparison_combined.png")

# ============================================================================
# PLOT 3: Improvement bar chart
# ============================================================================
plt.figure(figsize=(12, 7))

improvements = []
selection_labels = []

for selection_type in selection_types:
    best_scores = results_by_selection[selection_type]['best_scores']
    initial = best_scores[0]
    final = best_scores[-1]
    improvement = initial - final
    improvements.append(improvement)
    selection_labels.append(selection_type)

bars = plt.bar(range(len(selection_labels)), improvements, 
               color=colors, alpha=0.75, edgecolor='black', linewidth=1.5)
plt.xticks(range(len(selection_labels)), selection_labels, rotation=30, ha='right', fontsize=11)
plt.xlabel('Selection Type', fontsize=13, fontweight='bold')
plt.ylabel('Fitness Improvement', fontsize=13, fontweight='bold')
plt.title('Total Fitness Improvement by Selection Type', fontsize=15, fontweight='bold', pad=20)
plt.grid(True, alpha=0.3, axis='y', linestyle='--')

# Add value labels on bars
for i, bar in enumerate(bars):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selection_comparison_improvement.png'), dpi=300, bbox_inches='tight')
print("✓ Plot saved: selection_comparison_improvement.png")

# ============================================================================
# PLOT 3: Meilleur fitness comparaison (bar chart)
# ============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

final_best = [results_by_selection[st]['best_scores'][-1] for st in selection_types]
bars1 = ax1.bar(range(len(selection_labels)), final_best, 
                color=colors, alpha=0.75, edgecolor='black', linewidth=1.5)
ax1.set_xticks(range(len(selection_labels)))
ax1.set_xticklabels(selection_labels, rotation=30, ha='right', fontsize=10)
ax1.set_xlabel('Type de selection', fontsize=12, fontweight='bold')
ax1.set_ylabel('Score Fitness finale ', fontsize=12, fontweight='bold')
ax1.set_title(' Meilleur Fitness finale par type de Selection', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y', linestyle='--')

for i, bar in enumerate(bars1):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}',
            ha='center', va='bottom', fontsize=8, fontweight='bold')

final_worst = [results_by_selection[st]['worst_scores'][-1] for st in selection_types]
bars2 = ax2.bar(range(len(selection_labels)), final_worst, 
                color=colors, alpha=0.75, edgecolor='black', linewidth=1.5)
ax2.set_xticks(range(len(selection_labels)))
ax2.set_xticklabels(selection_labels, rotation=30, ha='right', fontsize=10)
ax2.set_xlabel('Type de selection', fontsize=12, fontweight='bold')
ax2.set_ylabel('Fitness finale Score', fontsize=12, fontweight='bold')
ax2.set_title('Pire Fitness finale par type de Selection', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y', linestyle='--')

for i, bar in enumerate(bars2):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}',
            ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selection_comparison_final.png'), dpi=300, bbox_inches='tight')
print("✓ Plot saved: selection_comparison_final.png")

# ============================================================================
# PRINT STATISTIQUES
# ============================================================================
print("\n" + "="*80)
print("STAT")
print("="*80)

for selection_type in selection_types:
    results = results_by_selection[selection_type]
    best_scores = results['best_scores']
    worst_scores = results['worst_scores']
    
    initial_best = best_scores[0]
    final_best = best_scores[-1]
    initial_worst = worst_scores[0]
    final_worst = worst_scores[-1]
    
    improvement_best = initial_best - final_best
    improvement_best_pct = (improvement_best / initial_best) * 100 if initial_best != 0 else 0
    
    improvement_worst = initial_worst - final_worst
    improvement_worst_pct = (improvement_worst / initial_worst) * 100 if initial_worst != 0 else 0
    
    print(f"\n{selection_type}:")
    print(f"  Meilleur individu:")
    print(f"    Initial:     {initial_best:.6f}")
    print(f"    Final:       {final_best:.6f}")
    print(f"  Pire Individu:")
    print(f"    Initial:     {initial_worst:.6f}")
    print(f"    Final:       {final_worst:.6f}")

# Trouver le meilleur type de selection 
best_selection = min(selection_types, key=lambda st: results_by_selection[st]['best_scores'][-1])
print("\n" + "="*80)
print(f"BEST SELECTION TYPE: {best_selection}")
print(f"Final fitness: {results_by_selection[best_selection]['best_scores'][-1]:.6f}")
print("="*80)

plt.show()
