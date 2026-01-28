# Algorithmes génétique pour l'Optimisation de Tables de Rotation d'ADN
https://github.com/Mc-cloud/Algo-genetique

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code Coverage](https://img.shields.io/badge/coverage-check%20tests-brightgreen.svg)](tests/)


## Objectif

Ce projet utilise un algorithme génétique pour optimiser les paramètres de rotation (twist, wedge, direction) des dinucléotides afin de minimiser la distance de fermeture des structures ADN circulaires (plasmides). L'algorithme cherche à trouver une table de rotations qui permette à la séquence ADN de se refermer sur elle-même avec une erreur minimale.

## Principe 

L'ADN est représenté par une trajectoire 3D calculée à partir : 
- D'une **séquence de nucléotides** (A, T, G, C)
- D'une **table de rotations** définissant trois angles pour chaque dinucléotide:
    - Twist (rotation autour de l'axe)
    - Wedge (inclinaison)
    - Direction (orientation)

L'algorithme génétique optimise cette table pour que la structure 3D forme un cercle fermé

## Fonctionnalités

- **Encodage génétique** : Table de rotations comme ADN de l'individu
- **Multiples méthodes de sélection** : 
    - Elitiste
    - Tournoi
    - Roulette
    - Par rang
- **Opérateurs génétiques** : 
    - Croisement pondéré par le fitness
    - Mutations (petites fréquentes et des plus gross rares)
    - Elitisme
- **Fonction de fitness** : Test de fermeture à plusieurs points de coupure
- **Recuit simulé** : Température décroissante pour affiner la convergence
- **Visualisations** : Génération de graphiques et GIfs d'évolution
- **Benchmarks** : Comparaison systématique de configurations


## Structure du projet

```
Algo-genetique/
├── algo/                      # Cœur de l'algorithme génétique
│   ├── algogenetique.py      # Classe Individu et fonction AlgoGenetique
│   ├── fitness.py            # Calcul du score de fermeture
│   └── selection.py          # 7 méthodes de sélection différentes
├── dna/                       # Représentation de l'ADN
│   ├── RotTable.py           # Table de rotations des dinucléotides
│   └── Traj3D.py             # Calcul de trajectoire 3D
├── data/                      # Séquences ADN de test
│   ├── plasmid_2k_*.fasta    # Plasmides de 2000 paires de bases
│   ├── plasmid_8k.fasta      # Plasmide de 8000 paires de bases
│   └── plasmid_180k.fasta    # Grand plasmide
├── data_algo/                 # Résultats d'expériences sauvegardés
├── documents/                 # Documentation et rapports
│   └── Rapport_*.pdf         # Rapport détaillé du projet
├── gifs/                      # Visualisations animées
│   ├── benchmark_*.gif       # Résultats de benchmarks
│   └── etapes.gif            # Évolution d'une simulation
├── tests/                     # Tests unitaires
│   ├── test_algogenetique.py
│   ├── test_fitness.py
│   └── test_selection.py
├── main.py                    # Script principal
├── plot.py                    # Génération de graphiques
├── resultsmanager.py          # Gestion des résultats
├── simulsmanager.py           # Gestion des simulations
├── executeur_comparaison_algos.py  # Comparaison d'algorithmes
├── benchmark.py               # Recherche de paramètres optimaux
├── benchmark_cuts.py          # Benchmark sur les points de coupure
└── tests_param.py             # Tests paramétriques
```

## 🔧 Installation
### En utilisant pip :
```bash
git clone https://github.com/Mc-cloud/Algo-genetique.git
cd Algo-genetique

pip install -r requirements.txt
```
### Avec Conda

```bash
# Cloner le dépôt
git clone https://github.com/Mc-cloud/Algo-genetique.git
cd Algo-genetique

# Créer et activer l'environnement conda
conda env create -f environment.yaml
conda activate algo-genetique
```

## Utilisation :
Mettre l'explication pour execution...

## Benchmarks 
```bash
# Benchmark avec recherche de grille automatique
python benchmark.py

# Tester l'impact du nombre de coupures
python benchmark_cuts.py

# Tests paramétriques personnalisés
python tests_param.py
```

## Utilisation

Le **fichier executeur** (`executeur_comparaison_algos.py`) est un programme interactif qui guide l'utilisateur pas à pas pour configurer et exécuter l'algorithme génétique.

```bash
python executeur_comparaison_algos.py
```

### Étape 1 : Sélection des fichiers d'entrée

Le programme vous demandera d'abord les fichiers nécessaires :

```
Indiquez le fichier '.fasta' contenant la séquence du plasmide d'étude.
> data/plasmid_8k.fasta
```

Fichiers FASTA disponibles dans `data/` :
- `plasmid_2k_*.fasta` : Petits plasmides (tests rapides)
- `plasmid_8k.fasta` : Plasmide de taille moyenne (recommandé)
- `plasmid_180k.fasta` : Grand plasmide (calculs longs)

```
Indiquez le fichier '.json' correspondant à la table de Rotation initiale.
S'il s'agit de la table du modèle, faites simplement 'Enter'
> [Enter pour utiliser dna/table.json par défaut]
```

### Étape 2 : Configuration des populations

Vous pouvez configurer **plusieurs populations** avec des paramètres différents pour les comparer :

```
Voulez-vous ajouter une population ?
Actuellement 0 populations prévues.
    oui/o
    non/n
> oui
```

Pour chaque population, vous devrez configurer :

#### a. Paramètres de fitness

```
Sur combien de bases voulez-vous tester la qualité du recollement ?
    • 1 ≤ n ≤ longueur(séquence ADN)
    • Par défaut n = 2
> 2
```
→ Nombre de nœuds à comparer entre début et fin (`nbappend`)

```
Combien d'autres points de départ voulez-vous tester ?
    • 0 ≤ n
    • Par défaut n = 0
> 0
```
→ Nombre de coupures supplémentaires (`nbcuts`). ⚠️ Attention : augmente le temps de calcul !

#### b. Méthode de sélection

```
Quelle façon de sélectionner les survivants ?
    1 élitiste
    2 tournoi
    3 roulette fitness
    4 roulette rang
    5 roulette rang géométrique
    6 roulette exponentielle
> 2
```

Ou tapez le nom complet : `tournoi`, `élitiste`, etc.

#### c. Taux de sélection

```
Quelle proportion de la population doit subsister ?
    • 0 < q < 1
    • Par défaut q = 0.5
> 0.3
```
→ Exemple : 0.3 signifie que 30% des individus deviennent géniteurs

#### d. Dimensionnement

```
Combien d'individus par génération ?
    • Par défaut n = 100
> 150
```

```
Combien de générations ?
    • Par défaut n = 20
> 50
```

### Étape 3 : Exécution et visualisation

Une fois toutes les populations configurées, l'exécuteur :

1. **Lance les simulations** séquentiellement
2. **Affiche la progression** en temps réel :
   ```
   Lancement de la 1-e population :
   itération : 1 / 50
   fit : 12.456
   Meilleur pour iter 1 : 3.234
   Pire pour iter 1 : 45.678
   ```

3. **Génère un graphique comparatif** montrant l'évolution du fitness de toutes les populations

## Tests :

Exécutez la suite de tests pour vérifier l'implémentation :

```bash
# Exécuter avec couverture
coverage run -m unittest discover -s tests -p "test_*.py"
coverage report
```
## Méthodes de sélection : 

Le projet implémente 7 méthodes de sélection :

| Méthode | Description | Usage |
|---------|-------------|-------|
| **elitiste** | Garde les N meilleurs individus | Convergence rapide, risque de convergence prématurée |
| **tournament** | Tournoi entre paires + 10% d'élite | Bon équilibre exploration/exploitation |
| **roulette** | Probabilité inversement proportionnelle au score | Maintient la diversité |
| **roulette_exp** | Roulette avec distribution exponentielle | Pression de sélection ajustable |
| **roulette_exp_norm** | Roulette exponentielle normalisée | Bon pour la convergence finale |
| **rang_reel** | Probabilité proportionnelle au rang | Évite la domination excessive |
| **rang_geo** | Distribution géométrique des probabilités | Bon compromis pression/diversité |

## Fonction de fitness

1. **Calcul de la trajectoire 3D**: Chaque dinucléotide applique une rotation
2. **Test de fermeture** : Calcule la distance euclidienne entre le début et la fin
3. **Multi-points** : teste à plusieurs points de coupure pour robustesse
4. **Score final** : Norme euclidienne des distances

## Visualisations

Les visualisations incluent : 
- Evolution du meilleur score
- Evolution du pire score
- Diversité de la population
- Comparaison entre méthodes

## Tests

```bash
# Exécuter tous les tests
python -m unittest discover -s tests

# Avec couverture de code
coverage run -m unittest discover -s tests -p "test_*.py"
coverage report

# Tests spécifiques
python -m unittest tests.test_fitness
python -m unittest tests.test_selection
python -m unittest tests.test_algogenetique
```

## Auteurs
- **Matheo Cahitte** [Mc-cloud](https://github.com/Mc-cloud)
- **Clément Cournil-Rabeux** 
- **Melkior Demaille**
- **Clément Rebola**

## Documentation

Pour plus de détails, consultez :
- Le rapport complet dans `documents/Rapport_EI_*.pdf`
- Les présentations dans `documents/AG-Pres.pdf` et `documents/AG-Poly.pdf`
- Le sujet initial dans `documents/EI_AlgoGen_project.pdf

**Statut** : Projet terminé (2024-2025)
