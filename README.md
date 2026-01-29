# Algorithmes génétique pour l'Optimisation de Tables de Rotation d'ADN
https://github.com/Mc-cloud/Algo-genetique

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code Coverage](https://img.shields.io/badge/coverage-check%20tests-brightgreen.svg)](tests/)


A cause de bugs sur gitlab, nous avons décidé de cloner le git vers github. Ainsi, la contribution de certains membres du groupes n'a pas été pris en compte par github car ils se sont connectés à github avec leur compte gitlab : Clement Cournil-Rabeux et Clement Rebola. 
Voici le lien vers le gitlab qui montre clairement qu'ils ont contribué énormément au projet : https://gitlab-cw2.centralesupelec.fr/clement.cournil-rabeux/jeux-evolutionnaires/
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
├── executeur_comparaison_algos.py  # Exécuteur d'algorithmes, 
interface dans terminal
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
→ Lors qu'une option `Par défaut` est proposée, il suffit de renvoyer un champ vide pour la sélectionner.

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

(Vous pouvez en utiliser d'autres, tant que vous précisez bien le `PATH`.)

Ensuite :
```
Indiquez le fichier '.json' correspondant à la table de Rotation initiale.
S'il s'agit de la table du modèle, faites simplement 'Enter'
> [Enter pour utiliser dna/table.json par défaut]
```
Laissé en option même si à priori la table de départ sera forcément celle du modèle.

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
→ Nombre de bases aux extrémités dont on teste la bonne superposition (`nbappend`).

```
Combien d'autres points de départ voulez-vous tester ?
    • 0 ≤ n
    • Par défaut n = 0
> 0
```
→ Nombre de coupures supplémentaires du plasmide testées (`nbcuts`), réparties de manière homogène. ⚠️ Attention : augmente le temps de calcul !

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

Indiquez le numéro ou tapez le nom complet : `tournoi`, `élitiste`, etc.

#### c. Taux de sélection

```
Quelle proportion de la population doit subsister ?
    • 0 < q < 1
    • Par défaut q = 0.3
> 0.4
```
→ Exemple : 0.4 signifie que 40% des individus sont conservés d'une génération à l'autre (parmis eux seront les géniteurs).

#### d. Dimensionnement

```
Combien d'individus par génération ?
    • Par défaut n = 100
> 150
```
Des résultats sensiblement meilleurs sont souvent obtenus avec une plus large population, par exemple `1000` individus.
```
Combien de générations ?
    • Par défaut n = 20
> 50
```
Évidemment le temps d'exécution est croissant avec ces deux paramètres.

### Étape 3 : Exécution, visualisation et téléchargement

#### a. Exécution

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
Attention, les seed d'aléatoire seront différentes pour chaque population, il peut être intéressant de comparer plusieurs instances aux paramètres identiques !

#### b. Visualisation de la fitness

Une fois les calculs terminés, l'exécuteur propose de **Générer des graphiques comparatifs** montrant l'évolution du fitness du meilleur de chaque populations, avec différents paramètres :

```
Voulez-vous afficher afficher l'évolution de la fitness du meilleur candidat de chaque population 
        • non/n 
        • selon un couple (n_bases_recollement,n_coupures) 
        • selon la fonction de fitness de chacune des population d'indice i_1,i_2,…,i_n ; i_k ≥ 1 
        • selon la fonction de fitness de chacune des populations choisies : t/tout 
>(2,1)
```
Ici il s'agit de choisir la (ou les) fonctions de fitness à employer pour comparer les meilleurs candidats des différentes populations.
- Indiquez directement un couple `(nbappend,nbcoup)` pour un seul graphique selon la fitness ayant ces paramètres.
- Indiquez une liste (par exemple `1, 3, 4`) des indices des populations dont on veut utiliser la fonction de fitness pour la comparaison. L'exemple donnera ici trois graphiques, où les meilleurs de *toutes* les populations seront comparées selon les *paramètres de fitness* de la première, la troisième et la quatrième population respectivement.
- Indiquez `t` ou `tout` pour afficher un graphique pour les paramètres de fitness de chaque population (cette option est donc équivalente au fait d'énumérer `1,2,3,4` s'il y a quatre populations au total, par exemple).

```
Appuyez sur Entrée pour fermer les graphiques et passer à l'étape suivante
```
Plusieurs figures peuvent apparaître, vous avez l'option de les télécharger directement avec l'interface de la fenêtre. Répondre à cet input les fermera.

#### c. Slider d'évolution des meilleures trajectoires

Le programme va demander si vous voulez afficher un widget de type "slider" pour visualiser dynamiquement l'évolution de la trajectoire du meilleur individu le long des génération. Attention, les sliders s'affichent l'un après l'autre.

```
Voulez-vous afficher l'évolution du meilleur chemin en fonction de la génération, pour chaque population ? 
        • o/oui 
        • n/non 
>o
```
→ Si vous acceptez, le programme affichera un widget slider *pour chaque population*, l'un après l'autre.
 
⚠️ Attention : Nous avons remarqué que Slider (issu de la bibliothèque matplotlib.widget) est parfois peu interactif, il est possible que sur certains ordinateurs il soit compliqué de faire glisser le curseur.

```
Appuyez sur Entrée pour fermer les graphiques et passer à l'étape suivante
```

De même, plusieurs figures vont apparaître. Malheureusement les télécharger ici ne sauvgardera qu'un arrêt sur image. Pour télécharger sous forme de `.gif` voir la fonction `save_trajectory_gif` de `plot.py`.

#### d. Téléchargement des tables de rotation des meilleurs individus.

Une fois toutes ces étapes passées, le programme va proposer de garder en mémoire les tables des meilleurs candidats des différentes populations, avec ces options :

```
Voulez-vous enregistrer les tables json du meilleur candidat : 
        • non/n 
        • seulement des populations d'indice i_1,i_2,…,i_n ; i_k ≥ 1 
        • de toutes les populations : t/tout 
>t
```
→ Comme pour la visualisation comparée des fitness de chaque population, on peut ici choisir une liste d'indices ( par exemple `1,3`) ou directement l'option `t`/`tout`.

Si au moins un indice est sélectionné :

```
Dans quel dossier enregistrer les json ? 
        (défaut: 'rot_tables_results') 
>
```
→ Le fichier par défaut est dédié à ces résultats. Vous pouvez en proposer un autre, tant qu'il existe et que son `PATH` est correcte.

Si tout est en ordre, le(s) fichier(s) sera/ont téléchargé(s) sous le format :
`[nom_dossier_résultat]/optimal_rot_table_[chemin_source_de_la_séquence]_[nom_du_fichier_de_la_séquence]_[méthode_de_sélection]_[nb_append]_[nb_cuts]_[nb_individus]_[nb_générations]`

Le programme affichant le chemin dans le terminal ; exemple :

```
rot_tables_results/optimal_rot_table_data_plasmid_8k.fasta_elitiste_4_3_150_25
```

### Fin
Le programme affiche
```
Programme terminé.
```

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

## Visualisations : 
Le projet génère diverses visualisations pour analyser les dynamiques évolutives :
- **Evolution du Fitness** : Suivi du fitness de la population au fil du temps
![png](evolutio_metrique.png)
- **Evolution d'un plasmide au fils des générations:**
![gif](gifs/ultimate.gif)
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
