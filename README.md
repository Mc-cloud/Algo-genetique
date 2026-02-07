# Genetic Algorithms for DNA Plasmid Optimization
https://github.com/Mc-cloud/Algo-genetique

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code Coverage](https://img.shields.io/badge/coverage-check%20tests-brightgreen.svg)](tests/)

**Genetic Algorithm** is a scientific computing tool designed to optimize the 3D folding of circular DNA structures. By adjusting the rotation parameters (twist, wedge, direction) of dinucleotides, the algorithm miniùizes the closure distance, ensuring the DNA forms a stable, closed loop.

## Project Overview

The project simulates DNA as a 3D trajectory based on : 
* **Nucleotide Sequence** (A, T, C, G)
* **Rotation Table:** Defines the geometric properties for every possible dinucleotide pair

## Key Features

* **Genetic Encoding:** Rotation tables are treated as the "DNA" of the solution.
* **Advanced Genetic Algorithm:** Implements 7 distinct selection strategies including **Elitist**, **Tournament**, **Roulette** (Linear/Exponential), and **Rank-based** selection.
* **3D Biological Modeling:** Uses a custom fitness function that calculates the **3D trajectory** of the DNA sequence and minimizes the **Euclidean closure distance** (gap between start and end).
* **Operators:** Weighted Crossover and Dynamic Mutation (simulated annealing).
* **Robustness:** Multi-objective optimization capabilities (optimizing closure at multiple internal cut points).
* **Interactive UI:** Built with **Streamlit** for real-time configuration and visualization.
* **Performance:** High-performance matrix calculations using `numba` JIT compilation.

## 📚 Documentation & Technical Report

For a deep dive into the mathematical models, fitness function definitions, and the convergence analysis of different selection strategies, please refer to the full technical report included in this repository:

* 📄 **[Technical Report (PDF)](documents/Rapport_EI_ST2_V2_ENG.pdf)**: Detailed explanation of the genetic algorithm implementation, parameter tuning, and performance benchmarks.

## Installation

### Using Conda
```bash
git clone [https://github.com/Mc-cloud/Algo-genetique.git](https://github.com/Mc-cloud/Algo-genetique.git)
cd Algo-genetique
conda env create -f environment.yaml
conda activate algo-genetique
```

### Using pip
```bash
pip install -r requirements.txt
```

## Usage
```bash
streamlit run scripts/app_client.py
```
* **Input Data:** Select a plasmid Fasta file
* **Population Config:** Define population size and generations
* **Selection Strategy:** Choose a strategy and selection rate
* **Number of cuts and added points:** IIf you want to minimize the distance, then you should add 0 cut. Adding cuts forces the algorithm to verify closure at multiple internal points (increasing robustness but computational cost). The added points (overlap) determine how many bases are repeated at the end to ensure the closure angle is smooth (tangent continuity); usually 1 or 2 is optimal.
* **Visualization:** Generate comparison plots and export the optimal rotation tables to json


## Benchmarks
You can run automated benchmarks to compare selection methods or parameter sensitivity

```bash
python scripts/benchmark.py

python scripts/benchmark_cuts.py
```

## Project Structure
```
Algo-genetique/
├── src/genetic_algo/      # Source code
│   ├── core/              # GA Logic (Individual, Selection, Fitness)
│   ├── dna/               # DNA 3D Modeling (RotTable, Traj3D)
│   └── utils/             # Visualization & I/O
├── data/                  # FASTA sequences
├── tests/                 # Unit tests
├── documents/             # Academic reports and presentations
└── scripts/               # Executable benchmarks and runners
```



## Visualisations : 
The project generates multiple visualizatins to analyse the dynamics of evolution : 
- **Evolution du Fitness** : Suivi du fitness de la population au fil du temps
![png](documents/images/evolutio_metrique.png)
- **Evolution d'un plasmide au fils des générations:**
![gif](documents/images/ultimate.gif)


## Tests
To execute the tests, please make sure to be in algo-genetique folder
```bash
# Execute all the tests
python -m unittest discover -s tests

# Shows coverage of the code
coverage run -m unittest discover -s tests -p "test_*.py"
coverage report

# Specific tests
python -m unittest tests.test_fitness
python -m unittest tests.test_selection
python -m unittest tests.test_algogenetique
```

## Autors
- **Matheo Cahitte** [Mc-cloud](https://github.com/Mc-cloud)
- **Clément Cournil-Rabeux** 
- **Melkior Demaille**
- **Clément Rebola**

**Status** : Project over (2026)
