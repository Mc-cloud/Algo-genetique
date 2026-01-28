"""
L'idée est d'avoir un programme interactif qui demande au client ce qu'il veut faire dans le terminal.
il doit pouvoir choisir s'il veut 
    - p0 choisir sa table de rotation et le .fasta à tester
    - p0 pouvoir faire un test
    - p1 choisir un seul mode de résolution et obtenir le ( p3 ou les k) meilleur(s) résultat(s)
    - p2 choisir plusieurs modes de résolution et obtenir leurs meilleurs résultats
    - p3 choisir ou non de récupérer les courbes d'évolution de la fitness, selon quelle fonction ( p4 et de représenter ou non les k meilleurs/pires/etc)
    - p2 choisir s'il veut évaluer les différents résultats, selon quelle fitenss
    - p1.5 permettre CTRL+Z pour revenir à la question précédente

"""
from json import load as json_load
from json import dump as json_dump
from dna.RotTable import *
from dna.Traj3D import *
from algo.algogenetique import AlgoGenetique
from algo.selection import selections_dic
from algo.fitness import fitness
from plot import *
import numpy as np
from simulsmanager import *
import time
import os

def checkpath(fichier):
    test = os.path.exists(fichier)
    if  not test:
            print(f"Erreur : Le fichier '{fichier}' est introuvable.")
            print("Vérifiez que :")
            print("  - Le fichier est dans le même dossier que ce programme")
            print("  - Ou indiquez le chemin complet (ex: /chemin/vers/fichier.txt)")
            print("  - Ou le chemin relatif (ex: ../dossier/fichier.txt)")
            print()
    return test

selections_dic ={
    "elitiste":"elitiste","tournament":"tournament","roulette":"roulette","rang_reel":"rang_reel","rang_geo":"rang_geo","roulette_exp":"roulette_exp",
    "élitiste":"elitiste","tournoi":"tournament", "roulette fitness":"roulette","roulette rang":"rang_reel","roulette rang géométrique":"rang_geo","roulette exponentielle":"roulette_exp",
    "1":"elitiste","2":"tournament","3":"roulette","4":"rang_reel","5":"rang_geo","6":"roulette_exp",
    "":"elitiste",
}


if __name__ == "__main__" :
    print("Logiciel de détermination de matrice de rotation optimale pour un plasmide")

    while True : # Demande le fichier lié à la séquence ADN
        sequence_file = input("Indiquez le fichier '.fasta' contenant la séquence du plasmide d'étude. \n> ")
        if not checkpath(sequence_file) :
            pass
        elif not sequence_file[-6:]==".fasta":
            print("Erreur : Le fichier doit avoir l'extension '.fasta'.")
            print()
        else:
             break
    print()
    # Read file
    lineList = [line.rstrip('\n') for line in open(sequence_file)]
    # Formatting
    seq = ''.join(lineList[1:])

    while True : # Demande le fichier lié à la table de rotation
            table_rot_file = input("Indiquez le fichier '.json' correspondant à la table de Rotation initiale sur laquelle itérer. \n S'il s'agit de la table du modèle, faites simplement 'Enter' \n>")
            if table_rot_file == "":
                table_rot_file = "dna/table.json"
                break
            if not checkpath(table_rot_file):
                pass
            elif not table_rot_file[-5:] == ".json":
                print("Erreur : Le fichier doit avoir l'extension '.json'.")
                print()
            else:
                break
    print()
    # while True :
    #Rajouter d'autres options comme décrit en haut de ce fichier (type de comparaison etc)

    list_nb_individus = []
    list_nb_gen = []
    list_taux_selec = []
    list_selection_type = []
    list_nb_append = []
    list_nb_cuts = []
    list_nb_individus = []
    list_poisson = []
    nb_instances = 0
    while True : # Demander la liste des algorithmes à exécuter, avec quels paramètres.
        new_gen = input(f"Voulez-vous ajouter une population (instance de l'algorithme dont il faudra préciser les paramètres) ? \nActuellement {nb_instances} populations prévues.  \n\toui/o \n\tnon/n \n>")
        if not new_gen in {'o','oui','n','non',''} :
            print("Erreur : veuillez donner une entrée valide.")
            print()
            pass
        elif new_gen in {'non','n',''}:
            break
        else:
            nb_instances+=1
        print()

        print("Choix des paramètres de la fonction de fitness.")
        while True : # paramètre de la fonction de fitness : nombre de recollements à tester
            var_nb_append = input("Sur combien de bases voulez-vous tester la qualité du recollement des extrémités du plasmide (qualité de la boucle) ? \n\t• 1 ≤ n ≤ longueur(séquence ADN). \n\t• Par défaut n = 2. \n Attention à ne pas dépasser la taille de la séquence ! \n>")
            if var_nb_append =="":
                var_nb_append = 2
                break
            elif not var_nb_append.isdigit() or not int(var_nb_append)>0:
                print("Entrez un nombre entier strictement positif en base décimale.")
                print()
            else:
                break
        print()
        list_nb_append.append(max(1,int(var_nb_append)))

        while True : # paramètre de la fonction de fitness : nombre de coupures supplémentaires
            var_nb_cuts = input("Combien d'autres points de départ du plasmide que celui donné par le fichier voulez-vous tester à chaque itération ? \n\t• 0 ≤ n \n\tPar défaut n = 0 \nAttention, le temps de calcul est proportionnel à 1 + ce facteur ! \n>")
            if var_nb_cuts =="":
                var_nb_cuts = 0
                break
            elif not var_nb_cuts.isdigit() or not int(var_nb_cuts)>=0:
                print("Entrez un nombre entier positif en base décimale.")
                print()
            else:
                break
        print()
        list_nb_cuts.append(max(1,int(var_nb_cuts)))

        print("Type de sélection.")
        while True : # Demander le type de selection
            selection_type_n = input("Quelle façon de sélectionner les survivants sur chaque génération ? \n\t1 élitiste \n\t2 tournoi \n\t3 roulette fitness\n\t4 roulette rang \n\t5 roulette rang géométrique \n\t6 roulette exponentielle \n>")
            if not selection_type_n in selections_dic:
                print("Erreur : entrez soit le chiffre correspondant à l'option, soit le nom de l'option.")
                print()
            else:
                break
        print()
        list_selection_type.append(selections_dic[selection_type_n])

        while True : # Demander la proportion des survivants
            var_taux_selec = input("Quel proportion de la population doit subsister d'une itération à l'autre ? \n\t• 0 < q < 1 \n\tPar défaut q = 0.5 \n>")
            if var_taux_selec == "":
                taux_selec = 0.5
                break
            try:
                taux_selec = float(var_taux_selec)
            except ValueError:
                print("Erreur : veuillez entrer un nombre décimal entre 0 et 1, le séparateur entre les parties entières et fractionnaires est le '.'.")
                print()
                pass
            if not 0 < taux_selec and taux_selec < 1 :
                print("Le taux de sélection doit être compris entre 0 et 1 strictement.")
                print()
            else:
                break
        print()
        list_taux_selec.append(taux_selec)

        print("Dimensionnement du processus.")
        while True : # dimension de la population initiale
            var_nb_indiv = input("Combien d'individus doit-il y avoir par génération ? \n\t• Nombre entier strictement positif \n\t• par défaut n = 100 \n>")
            if var_nb_indiv =="":
                var_nb_indiv = 100
                break
            elif not var_nb_indiv.isdigit() or not int(var_nb_indiv)>0:
                print("Entrez un nombre entier strictement positif en base décimale.")
                print()
            else:
                break
        print()
        list_nb_individus.append(int(var_nb_indiv))

        while True : # nombre de générations
            var_nb_gen = input("Combien de générations doit-il y avoir pour cette exécution ? \n\t• Nombre entier strictement positif \n\t• par défaut n = 20 \n>")
            if var_nb_gen =="":
                var_nb_gen = 20
                break
            elif not var_nb_gen.isdigit() or not int(var_nb_gen)>0:
                print("Entrez un nombre entier strictement positif en base décimale.")
                print()
            else:
                break
        print()
        list_nb_gen.append(int(var_nb_gen))    
    # while True :
    #Rajouter d'autres options comme décrit en haut de ce fichier (autre endroit pour demander le type de comparaisons)
    print(f"Nombre total de populations : {nb_instances}")
    print()
    Best_indiv_list,Best_indiv_score_list,worst_indiv_score_list = [],[],[]
    for i in range(nb_instances):
        print(f"Lancement de la {i+1}-e population :")
        Best_indiv,Best_indiv_score,worst_indiv_score = AlgoGenetique(table_rot_file ,seq , nb_individus = list_nb_individus[i],nb_generations = list_nb_gen[i],taux_selec = list_taux_selec[i], selection_type = list_selection_type[i],poisson=False,nb_cuts = list_nb_cuts[i],nb_append = list_nb_append[i],recuit=False)
        print()
        Best_indiv_list.append(Best_indiv)
        Best_indiv_score_list.append([fitness(indiv.Rot_table,seq,nbcuts = 0) for indiv in Best_indiv])
        worst_indiv_score_list.append(worst_indiv_score)
        plt.plot(Best_indiv_score, label=f"{i+1}-ème population")
    plt.legend()
    plt.show()
    