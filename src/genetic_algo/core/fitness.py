import numpy as np
from dna.RotTable import RotTable
from dna.Traj3D import Traj3D


def dist_df(coords: list, nbappend = 1):
    """Prend en entrée une trajectoire sous la forme d'une liste d'array 
    correspondant aux coordonnées des nœuds, en joiniant les nbappend premiers
    nœuds à la fin. On compare ensuite la distance euclidienne entre les
    nœuds initiaux et ceux rajoutés."""

    """
    Calcule la distance de fermeture d'une trajectoire ADN.
    Mesure la distance euclidienne entre les premiers et derniers nœuds
    d'une trajectoire pour évaluer si la structure ADN se referme correctement.
    Args:
        coords: Liste d'arrays numpy représentant les coordonnées 3D des nœuds de la trajectoire
        nbappend: Nombre de nœuds à comparer entre le début et la fin (défaut: 1)
    Returns:
        Distance euclidienne totale entre les nbappend premiers et derniers nœuds
    Note:
        Utilisé pour vérifier la circularité de l'ADN (plasmides circulaires)
    """
    start = coords[:nbappend]
    end = coords[-nbappend:]
    distsq = 0

    for i in range(nbappend):
        distsq+= np.linalg.norm(start[i]-end[i])**2

    return np.sqrt(distsq)

def dist_euclid(scores: list) -> float:
    """
    Calcule la norme euclidienne d'un vecteur de scores.
    Args:
        scores: Liste de valeurs numériques
    Returns:
        Norme L2 (euclidienne) du vecteur de scores
    Note:
        Utilisé comme fonction de combinaison pour agréger plusieurs scores
    """
    return np.linalg.norm(scores)


def fitness(rot_table: RotTable, seq: str, fct_poids = dist_df, nbappend = 2, nbcuts = 2, coup_combin = dist_euclid) -> float :
    """
    Fonction de fitness pour évaluer la qualité d'une séquence ADN.
    
    Calcule un score de fermeture en testant la séquence à différents points de coupure
    pour évaluer sa capacité à former une structure circulaire stable.
    
    Args:
        rot_table: Table des rotations pour calculer la trajectoire 3D de l'ADN
        seq: Séquence ADN à évaluer (chaîne de caractères)
        fct_poids: Fonction de calcul du poids/score pour une trajectoire (défaut: dist_df)
        nbappend: Nombre de nœuds à ajouter à la fin pour tester la fermeture (défaut: 2)
        nbcuts: Nombre de points de coupure à tester (défaut: 2)
        coup_combin: Fonction pour combiner les scores des différentes coupures (défaut: dist_euclid)
    
    Returns:
        Score de fitness (plus le score est bas, meilleure est la fermeture)
    
    Note:
        La méthode teste plusieurs rotations circulaires de la séquence pour évaluer
        la stabilité globale de la structure ADN circulaire.
    """
    nbases = len(seq)
    assert nbases >= nbappend
    traj = Traj3D()

    def eval_une_coupure(seq: str, nbappend: int, indcut: int):
        """
        Évalue le score de fermeture pour un point de coupure donné.
        Args:
            seq: Séquence ADN
            nbappend: Nombre de nœuds à rajouter
            indcut: Index du point de coupure
        
        Returns:
            Score de fermeture pour cette coupure
        """
        traj.compute(seq[indcut:]+seq[:indcut+nbappend], rot_table)
        coords = traj.getTraj()
        score = fct_poids(coords,nbappend)
        return score
    
    list_coupes = [i*int(np.floor(nbases/(nbcuts+1))) for i in range(nbcuts+1)]

    score = coup_combin([eval_une_coupure(seq, nbappend, index) for index in list_coupes])

    return score

def fitness_basic(rot_table:RotTable, seq: str):
    return fitness(rot_table,seq,nbcuts=0)

#Tests

str_data = 'AACTGTCAGCTACCGATCATCTAGCTCTATATCGCGCATTAGCAGC'
rot_table = RotTable({
    "AA": [35.62 , 7.2 , -154 ,      0.06 ,  0.6   , 0],
    "AC": [34.4  , 1.1 ,  143 ,      1.3  ,  5     , 0],
    "AG": [27.7  , 8.4 ,    2 ,      1.5  ,  3     , 0],
    "AT": [31.5  , 2.6 ,    0 ,      1.1  ,  2     , 0],
    "CA": [34.5  , 3.5 ,  -64 ,      0.9  , 34     , 0],
    "CC": [33.67 , 2.1 ,  -57 ,      0.07 ,  2.1   , 0],
    "CG": [29.8  , 6.7 ,    0 ,      1.1  ,  1.5   , 0],
    "CT": [27.7  , 8.4 ,   -2 ,      1.5  ,  3     , 0],
    "GA": [36.9  , 5.3 ,  120 ,      0.9  ,  6     , 0],
    "GC": [40    , 5   ,  180 ,      1.2  ,  1.275 , 0],
    "GG": [33.67 , 2.1 ,   57 ,      0.07 ,  2.1   , 0],
    "GT": [34.4  , 1.1 , -143 ,      1.3  ,  5     , 0],
    "TA": [36    , 0.9 ,    0 ,      1.1  ,  2     , 0],
    "TC": [36.9  , 5.3 , -120 ,      0.9  ,  6     , 0],
    "TG": [34.5  , 3.5 ,   64 ,      0.9  , 34     , 0],
    "TT": [35.62 , 7.2 ,  154 ,      0.06 ,  0.6   , 0]
})

if __name__ == "__main__" :
    test1 = fitness(rot_table, str_data, nbappend= 3, nbcuts=2)
    test2 = fitness(rot_table, str_data, nbappend= 5, nbcuts=0)
    test3 = fitness(rot_table, str_data, nbappend= 3, nbcuts=5)
    print(test1)
    print(test2)
    print(test3)


    #test4 = fitness(rot_table, "AA", nbappend = 6, nbcuts=0)

    test5 = fitness(rot_table, "ATGC", nbappend = 2, nbcuts = 6)
    test6 = fitness(rot_table, "ATGC", nbappend = 2, nbcuts = 12)
    print(test5)
    print(test6)


    lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
    seq = ''.join(lineList[1:])
    test7 = fitness(rot_table,seq, nbappend = 2, nbcuts = 20)
    print(test7)