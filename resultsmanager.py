import pickle
import os

def save_simulation_data(filename, listes_data, params):

    payload = {
        "parameters": params,
        "data": {
            "liste_1": listes_data[0],
            "liste_2": listes_data[1],
            "liste_3": listes_data[2]
        }
    }
    
    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)

    try:
        with open(filename, 'wb') as f:
            pickle.dump(payload, f)
        print(f"Sauvegarde réussie dans '{filename}' avec les paramètres : {params}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")

def load_simulation_data(filename, check_nb_indiv, check_nb_gen, check_taux, check_type):
    """
    Charge les données SEULEMENT si les paramètres correspondent.
    
    Returns:
        tuple: (liste1, liste2, liste3) si tout est bon.
    Raises:
        ValueError: Si les paramètres ne correspondent pas.
        FileNotFoundError: Si le fichier n'existe pas.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Le fichier '{filename}' n'existe pas.")

    with open(filename, 'rb') as f:
        loaded_payload = pickle.load(f)

    stored_params = loaded_payload["parameters"]
    
    # 2. VÉRIFICATION DES PARAMÈTRES (Le "Check")
    errors = []
    
    if stored_params['nb_indiv'] != check_nb_indiv:
        errors.append(f"nb_indiv attendu {check_nb_indiv}, trouvé {stored_params['nb_indiv']}")
        
    if stored_params['nb_generations'] != check_nb_gen:
        errors.append(f"nb_generations attendu {check_nb_gen}, trouvé {stored_params['nb_generations']}")
        
    if stored_params['taux_selec'] != check_taux:
        errors.append(f"taux_selec attendu {check_taux}, trouvé {stored_params['taux_selec']}")
        
    if stored_params['selection_type'] != check_type:
        errors.append(f"selection_type attendu '{check_type}', trouvé '{stored_params['selection_type']}'")

    if errors:
        error_msg = "\n".join(errors)
        raise ValueError(f"❌ CONFLIT DE PARAMÈTRES :\n{error_msg}\nImpossible de charger des données incohérentes.")

    print("Paramètres validés. Chargement des données...")
    
    # 3. Retour des listes
    data = loaded_payload["data"]
    return data["liste_1"], data["liste_2"], data["liste_3"]