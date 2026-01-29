"""
L'idée est d'avoir un programme interactif qui demande au client ce qu'il veut faire dans le terminal.
il doit pouvoir choisir s'il veut 
    - p0 choisir sa table de rotation et le .fasta à tester
    - p0 pouvoir faire un test
    - p1 choisir un seul mode de résolution et obtenir le ( p3 ou les k) meilleur(s) résultat(s)
    - p2 choisir plusieurs modes de résolution et obtenir leurs meilleurs résultats
    - p3 choisir ou non de récupérer les courbes d'évolution de la fitness, selon quelle fonction ( p4 et de représenter ou non les k meilleurs/pires/etc)
    - p2 choisir s'il veut évaluer les différents résultats, selon quelle fitenss
    - p1 pouvoir enregistrer les .json des tables des meilleurs
    - p5 permettre CTRL+Z pour revenir à la question précédente

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
            print(f"Erreur : Le fichier/dossier '{fichier}' est introuvable.")
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
    print("Programme de détermination de matrice de rotation optimale pour un plasmide.")
    print()

    while True : # Demande le fichier lié à la séquence ADN
        sequence_file = input("Indiquez le fichier '.fasta' contenant la séquence du plasmide d'étude. \n> ")
        if not checkpath(sequence_file) :
            continue
        elif not sequence_file[-6:]==".fasta":
            print("Erreur : Le fichier doit avoir l'extension '.fasta'.")
            print()
            continue
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
                continue
            elif not table_rot_file[-5:] == ".json":
                print("Erreur : Le fichier doit avoir l'extension '.json'.")
                print()
                continue
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
        if not new_gen in {'o','oui','n','non','','y','yes','n','no'} :
            print("Erreur : veuillez donner une entrée valide.")
            print()
            continue
        elif new_gen in {'non','no','n',''}:
            break
        else:
            nb_instances+=1
        print()

        print("Choix des paramètres de la fonction de fitness.")
        print()
        while True : # paramètre de la fonction de fitness : nombre de recollements à tester
            var_nb_append = input("Sur combien de bases voulez-vous tester la qualité du recollement des extrémités du plasmide (qualité de la boucle) ? \n\t• 1 ≤ n ≤ longueur(séquence ADN). \n\t• Par défaut n = 2. \n Attention à ne pas dépasser la taille de la séquence ! \n>")
            if var_nb_append =="":
                var_nb_append = 2
                break
            elif not var_nb_append.isdigit() or not int(var_nb_append)>0:
                print("Entrez un nombre entier strictement positif en base décimale.")
                print()
                continue
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
                continue
            else:
                break
        print()
        list_nb_cuts.append(max(0,int(var_nb_cuts)))

        print("Type de sélection.")
        print()
        while True : # Demander le type de selection
            selection_type_n = input("Quelle façon de sélectionner les survivants sur chaque génération ? \n\t1 élitiste \n\t2 tournoi \n\t3 roulette fitness\n\t4 roulette rang \n\t5 roulette rang géométrique \n\t6 roulette exponentielle \n>")
            if not selection_type_n in selections_dic:
                print("Erreur : entrez soit le chiffre correspondant à l'option, soit le nom de l'option.")
                print()
                continue
            else:
                break
        print()
        list_selection_type.append(selections_dic[selection_type_n])

        while True : # Demander la proportion des survivants
            var_taux_selec = input("Quel proportion de la population doit subsister d'une itération à l'autre ? \n\t• 0 < q < 1 \n\tPar défaut q = 0.3 \n>")
            if var_taux_selec == "":
                taux_selec = 0.3
                break
            try:
                taux_selec = float(var_taux_selec)
            except ValueError:
                print("Erreur : veuillez entrer un nombre décimal entre 0 et 1, le séparateur entre les parties entières et fractionnaires est le '.'.")
                print()
                continue
            if not 0 < taux_selec and taux_selec < 1 :
                print("Le taux de sélection doit être compris entre 0 et 1 strictement.")
                print()
                continue
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
                continue
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
                continue
            else:
                break
        print()
        list_nb_gen.append(int(var_nb_gen))
    
    plot_bs = []
    print("Début des calculs")
    print(f"Nombre total de populations : {nb_instances}")
    print()
    Best_indiv_list = []

    for i in range(nb_instances):
        print(f"Lancement de la {i+1}-e population :")
        Best_indiv,Best_indiv_score,worst_indiv_score = AlgoGenetique(table_rot_file ,seq , nb_individus = list_nb_individus[i],nb_generations = list_nb_gen[i],taux_selec = list_taux_selec[i], selection_type = list_selection_type[i],poisson=False,nb_cuts = list_nb_cuts[i],nb_append = list_nb_append[i],recuit=False)
        print()
        Best_indiv_list.append(Best_indiv)
    print("Fin.")
    print()


    print("Options d'affichage")
    while True : # Options d'affichage : voulez-vous rajouter un affichage des fitness au cours du temps ?
        var_plot_bs = input("Voulez-vous afficher afficher l'évolution de la fitness du meilleur candidat de chaque population \n\t• non/n \n\t• selon un couple (n_bases_recollement,n_coupures) \n\t• selon la fonction de fitness de chacune des population d'indice i_1,i_2,…,i_n ; i_k ≥ 1 \n\t• selon la fonction de fitness de chacune des populations choisies : t/tout \n>")
        if var_plot_bs in {'','n','no','non'}:
            plot_bs = []
            break
        elif var_plot_bs in {'t','tout','a','all'}:
            plot_bs = range(nb_instances)
            break
        elif var_plot_bs[0]=='(':
            try:
                # Enlève parenthèses et remplace virgules par espaces
                fx_plot_bs = var_plot_bs.strip().strip('()').replace(',', ' ')
                fx_plot_bs = [int(x) for x in fx_plot_bs.split()]
                
                if len(fx_plot_bs) != 2:
                    print("Erreur : entrez exactement 2 nombres.\n")
                    continue
                else:
                    list_nb_append.append(fx_plot_bs[0])
                    list_nb_cuts.append(fx_plot_bs[1])
                    plot_bs = [nb_instances]
                    break 
            except ValueError:
                print("Erreur : format invalide. Exemple : (5, 10). Ne pas mettre de parenthèses pour une liste d'indices.\n")
                print()
                continue

        else:
            try:
                # Remplace virgules par espaces
                var_plot_bs = var_plot_bs.replace(',', ' ')
                # Convertit en liste d'entiers
                plot_bs = [int(x)-1 for x in var_plot_bs.split()]
                var_bool = True
                for i in plot_bs:
                    if i < 0 or i >= nb_instances:
                        print("N'indiquez que des numéros valides de population.")
                        var_bool = False
                        break
                if not var_bool:
                    continue
                else:
                    break
            except ValueError:
                print("Erreur : entrez une option correcte, 'n'/'non', 't/tout', un tuple parenthésé (n,k), ou des nombres entiers positifs séparés par des espaces ou des virgules, correspondant aux indices (i≥1) à choisir.")
                print()
                continue
    print()
    list_fitness_to_plot = {i : (lambda rot_table_x, idx=i: fitness(rot_table_x, seq, nbappend=list_nb_append[idx], nbcuts=list_nb_cuts[idx])) for i in plot_bs}
 

    for k in plot_bs:
        Best_indiv_score_list = []
        plt.figure(k+1)
        plt.title(f"Évolution de l'évaluation des meilleurs de chaque population pour la fonction de fitness {k+1}")

        for i in range(nb_instances):
            Best_indiv_score_list.append([list_fitness_to_plot[k](indiv.Rot_table) for indiv in Best_indiv_list[i]])
            plt.plot(Best_indiv_score_list[i], label=f"{i+1}-ème population")

        plt.legend()
        plt.show(block=False)
    
    if plot_bs:  # Si des plots ont été créés
        input("\nAppuyez sur Entrée pour fermer les graphiques et passer à l'étape suivante")
        print()
        plt.close('all')
    
    while True : # Demander s'il faut afficher les sliders pour les meilleurs résultats
        var_sliders = input("Voulez-vous afficher l'évolution de la trajectoire du meilleur de chaque génération avec un slider, pour chaque population ? \n(Parfois, le slider du widget est peu interactif.) \n\t• o/oui \n\t• n/non \n>")
        if not var_sliders in {'o','oui','n','non','','y','yes','n','no'} :
            print("Erreur : veuillez donner une entrée valide.")
            print()
            continue
        elif var_sliders in {'non','no','n',''}:
            break
        else:
            for indiv in Best_indiv_list:
                plot_with_slider(get_trajectories(indiv,seq), block=False)
            input("\nAppuyez sur Entrée pour fermer les graphiques et passer à l'étape suivante")
            print()
            break

    while True : # Demander s'il faut enregistrer les tables json, et pour quelle(s) population(s).
        var_download_json = input("Voulez-vous enregistrer les tables json du meilleur candidat : \n\t• non/n \n\t• seulement des populations d'indice i_1,i_2,…,i_n ; i_k ≥ 1 \n\t• de toutes les populations : t/tout \n>")
        if var_download_json in {'','n','no','non'}:
            download_json=[]
            break
        elif var_download_json in {'t','tout','a','all'}:
            download_json = range(nb_instances)
            break
        else:
            try:
                # Remplace virgules par espaces
                var_download_json = var_download_json.replace(',', ' ')
                # Convertit en liste d'entiers
                download_json = [int(x)-1 for x in var_download_json.split()]
                var_bool = True
                for i in download_json:
                    if i < 0 or i >= nb_instances:
                        print("N'indiquez que des numéros valides de population.")
                        var_bool = False
                        break
                if not var_bool:
                    continue
                else:
                    break
            except ValueError:
                print("Erreur : entrez une option correcte, 'n'/'non', 't/tout', ou des nombres entiers positifs séparés par des espaces ou des virgules, correspondant aux indices (i≥1) à choisir.")
                print()
                continue
       
    print()
    list_json_to_download = {i : f"optimal_rot_table_{sequence_file.replace('/', '_').replace('\\', '_')}_{list_selection_type[i]}_{list_nb_append[i]}_{list_nb_cuts[i]}_{list_nb_individus[i]}_{list_nb_gen[i]}"  for i in download_json}  

    while True : # Demander où enregistrer tout ça
        folder = input("Dans quel dossier enregistrer les json ? \n\t(défaut: 'rot_tables_results') \n>")
        if folder == "":
            folder = "rot_tables_results"
            break
        elif not checkpath(folder):
            continue
        else:
            break
    for i in list_json_to_download:
        filepath = os.path.join(folder, list_json_to_download[i])
        print(filepath)
        Best_indiv_list[i][-1].Rot_table.save(filepath)
        
        
    
    print("\nProgramme terminé.")