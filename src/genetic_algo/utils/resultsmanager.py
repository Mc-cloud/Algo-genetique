import pickle
import os

def save_simulation_data(filename, best_indiv_list, best_score_list, worst_score_list, params, dna_seq):
    """
    Sauvegarde 3 listes passées séparément.
    
    Args:
        filename (str): Chemin du fichier.
        best_indiv_list, best_score_list, worst_score_list : Tes trois listes de données.
        params (dict): Tes paramètres (nb_indiv, etc.).
        dna_seq (str): La séquence d'ADN.
    """

    full_params = params.copy()
    full_params['dna_seq'] = dna_seq

    payload = {
        "parameters": full_params,
        "data": {
            "best_indiv_list": best_indiv_list,
            "best_score_list": best_score_list,
            "worst_score_list": worst_score_list
        }
    }

    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)

    try:
        with open(filename, 'wb') as f:
            pickle.dump(payload, f)
        print(f"Sauvegarde réussie dans '{filename}'.")
    except Exception as e:
        print(f"Erreur sauvegarde : {e}")

def load_simulation_data(filename, check_dna_seq):
    """
    Charge les données , uniquement si l'adn est bien correspondant.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Fichier '{filename}' introuvable.")

    with open(filename, 'rb') as f:
        loaded_payload = pickle.load(f)

    stored_params = loaded_payload["parameters"]
    errors = []

    # --- VERIFICATIONS ---

    stored_dna = stored_params.get('dna_seq', "") 
    if stored_dna != check_dna_seq:
        disp_store = (stored_dna[:20] + '...') if len(stored_dna) > 20 else stored_dna
        disp_check = (check_dna_seq[:20] + '...') if len(check_dna_seq) > 20 else check_dna_seq
        errors.append(f"- ADN DIFFERENT ! Stocké: '{disp_store}' vs Demandé: '{disp_check}'")

    if errors:
        raise ValueError("PARAMÈTRES INCORRECTS :\n" + "\n".join(errors))

    print("Paramètres et ADN validés. Chargement...")
    d = loaded_payload["data"]
    return d["best_indiv_list"], d["best_score_list"], d["worst_score_list"],loaded_payload["parameters"]