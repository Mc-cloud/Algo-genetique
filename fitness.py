import numpy as np
import numpy.typing as npt
from dna.RotTable import RotTable
from dna.Traj3D import Traj3D


def dist_df(coords: list, nbappend = 1):
    """Prend en entrée une trajectoire sous la forme d'une liste d'array 
    correspondant aux coordonnées des nœuds, en joiniant les nbappend premiers
    nœuds à la fin. On compare ensuite la distance euclidienne entre les
    nœuds initiaux et ceux rajoutés."""
    start = coords[:nbappend]
    end = coords[-nbappend:]
    distsq = 0

    for i in range(nbappend):
        distsq+= np.linalg.norm(start[i]-end[i])**2

    return np.sqrt(distsq)



def fitness(rot_table: RotTable, seq: str, fct_poids = dist_df, nbappend = 1, nbcuts = 0) -> float :
    """prend en entrée une table de rotations et un code ADN,
    calcule le chemin et en déduit un score.
    Edit : prend aussi un nombre de nœuds finaux à rajouter, 
    et un nombre de coupures à faire."""
    nbases = len(seq)
    assert nbases >= nbappend
    
    def eval_une_coupure(seq: str, nbappend: int, indcut: int):
        traj = Traj3D()
        traj.compute(seq[indcut:]+seq[:indcut+nbappend], rot_table)
        coords = traj.getTraj()
        score = fct_poids(coords,nbappend)
        return score
    
    list_coupes = [i*np.floor(nbases/(nbcuts+1)) for i in range(nbcuts+1)]

    score_sq = sum(eval_une_coupure(seq, nbappend, index)**2 for index in list_coupes)

    return np.sqrt(score_sq)





