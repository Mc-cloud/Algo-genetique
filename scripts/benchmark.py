from src.genetic_algo.dna.RotTable import *
from src.genetic_algo.dna.Traj3D import *
from src.genetic_algo.core.algogenetique import AlgoGenetique
from src.genetic_algo.core.selection import selections_dic
from src.genetic_algo.core.fitness import fitness
import numpy as np
from src.genetic_algo.utils.simulsmanager import simul_and_save_results
from src.genetic_algo.utils.resultsmanager import load_simulation_data
from src.genetic_algo.utils.plot import plot_with_slider,get_trajectories,plot_best_worst,plot_three_indicators,save_trajectory_gif

def get_indicators(coords):
    dist = np.linalg.norm(coords[0]-coords[-2])
    v1 = coords[-1]-coords[-2]
    v2 = coords[1] -coords[0]
    v1,v2 = v1/np.linalg.norm(v1),v2/np.linalg.norm(v2)
    return dist,np.linalg.norm(v1-v2),np.dot(v2,v1)

'''
score = fitness(best.Rot_table,base_seq,nbcuts=0)
print(" score : ",best.score,"score final : ",score," via type de selection : ","elitiste")
traj_res = Traj3D(want_to_plot=True)
traj_res.compute(base_seq,best[-1].Rot_table)
traj_res.draw() #'''

def grid_search_params_save(base_save_filename,dna_seq,params_listed):
    """
    Docstring for grid_search_params
    
    :param params_listed: Dictionnaire avec en clef les paramètres usuels, et en valeur la liste de params à tester
    """
    if "nb_individus" not in params_listed:
        params_listed["nb_individus"]=150
    if "nb_generations" not in params_listed:
        params_listed[""]

    for nb_individus in params_listed.get("nb_individus",[150]):
        for nb_generations in params_listed.get("nb_generations",[20]):
            for taux_selec in params_listed.get("taux_selec",[0.5]):
                for selection_type in params_listed.get("selection_type",["elitiste"]):
                    for poisson in params_listed.get("poisson",[False]):
                        for nb_cuts in params_listed.get("nb_cuts",[0]):
                            for nb_append in params_listed.get("nb_append",[1]):
                                for recuit in params_listed.get("recuit",[False]):
                                    curr_params = {
                                        "nb_individus":nb_individus,
                                        "nb_generations":nb_generations,
                                        "taux_selec":taux_selec,
                                        "selection_type":selection_type,
                                        "poisson":poisson,
                                        "nb_cuts":nb_cuts,
                                        "nb_append":nb_append,
                                        "recuit":recuit
                                    }
                                    simul_and_save_results(base_save_filename+"_".join([f"{a}{curr_params[a]}" for a in curr_params]),dna_seq,curr_params)

def grid_search_compare(base_save_filename,dna_seq,params_listed,show="convergence_best",gif_name=""):
    """
    Docstring for grid_search_compare
    
    :param base_save_filename: chemin de base pour retrouver les simulations sauvegardées (préfixe commun à toutes les simulations)
    :param dna_seq: la séquence d'ADN à considérer
    :param params_listed: Dictionnaire avec en clef les paramètres usuels, et en valeur la liste de params à tester
    :param show: par défaut, show="convergence_best", affiche la convergence de la meilleure configuration. 
                 Si show="3indicators", affiche les 3 indicateurs de la meilleure configuration, étant l'alignement, la continuité et la distance.
                 sinon affiche la trajectoire de la meilleure configuration , en fonction des itérations.
    :param gif_name: vide par défaut, s'il n'est pas vide, sauvegarde un gif de l'évolution de la meilleure configuration sous ce nom.
    """
    best_config_score = float('inf')
    best_config=None
    best_res=None
    for nb_individus in params_listed.get("nb_individus",[150]):
        for nb_generations in params_listed.get("nb_generations",[20]):
            for taux_selec in params_listed.get("taux_selec",[0.5]):
                for selection_type in params_listed.get("selection_type",["elitiste"]):
                    for poisson in params_listed.get("poisson",[False]):
                        for nb_cuts in params_listed.get("nb_cuts",[0]):
                            for nb_append in params_listed.get("nb_append",[1]):
                                for recuit in params_listed.get("recuit",[False]):
                                    curr_params = {
                                        "nb_individus":nb_individus,
                                        "nb_generations":nb_generations,
                                        "taux_selec":taux_selec,
                                        "selection_type":selection_type,
                                        "poisson":poisson,
                                        "nb_cuts":nb_cuts,
                                        "nb_append":nb_append,
                                        "recuit":recuit
                                    }
                                    loc_path = base_save_filename+"_".join([f"{a}{curr_params[a]}" for a in curr_params])
                                    res = load_simulation_data(loc_path,dna_seq)
                                    indiv_list,b_list,w_list,_ = res
                                    final_score = b_list[-1]
                                    if final_score < best_config_score:
                                        best_config_score = final_score
                                        best_config = curr_params
                                        best_res = (indiv_list,b_list,w_list)
    print(f"La meilleur configuration trouvée est : {best_config}, avec un score de {best_config_score}")
    best_indiv_list,best_b_list,best_w_list = best_res
    def get_where_belong(param):
        return f"{best_config[param]}∈"+"{"+((f"{params_listed[param]}")[1:])[:-1]+"}" if len(params_listed[param])>1 else f"{best_config[param]}"
    title = f"\n nombre d'individus : {get_where_belong("nb_individus")} nombre d'itérations : {get_where_belong("nb_generations")} taux de sélection : {get_where_belong("taux_selec")} \nsélection:{get_where_belong("selection_type")} nombre de coupes : {get_where_belong("nb_cuts")} Poisson utilisé : {"Oui" if best_config["poisson"] else "Non"}"
    if show=="convergence_best":
        print("Visualisation de la convergence de la meilleure configuration : ")
        plot_best_worst(base_save_filename+"_".join([f"{a}{best_config[a]}" for a in best_config]),dna_seq,title="Evolution du pire et meilleur individu en fonction de l'itération, sur la meilleure configuration préalable de paramètres trouvée \n"+title)
    elif show=="3indicators":
        print("Visualisation des 3 indicateurs")
        dist,norm,ps = [],[],[]
        for indiv in indiv_list:
            traj_res = Traj3D()
            traj_res.compute(dna_seq,indiv.Rot_table)
            d,n,p = get_indicators(traj_res.getTraj())
            dist.append(d)
            norm.append(n)
            ps.append(p)
        plot_three_indicators(dist,norm,ps,title=title)
    else:
        print("Visualisation de la meilleure configuration..")
        plot_with_slider(get_trajectories(best_indiv_list,dna_seq))
    if gif_name:
        save_trajectory_gif(get_trajectories(best_indiv_list,dna_seq),gif_name,fps=10)

                                    

