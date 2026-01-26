import numpy as np
import numpy.typing as npt
from dna.RotTable import RotTable
from dna.Traj3D import Traj3D


def dist_df(coords: list):
    """Prend en entrée une trajectoire sous la forme d'une liste d'array 
    correspondant aux coordonnées des nœuds, garde uniquement premier et dernier
    pour évaluer leur distance euclidienne"""
    start, finish = coords[0], coords[-1]
    score = np.linalg.norm(finish-start)
    return score



def fitness(rot_table: RotTable, seq: str, fct_poids = dist_df, appendend = False) -> float :
    """prend en entrée une table de rotations et un code ADN,
    calcule le chemin et en déduit un score."""
    assert seq != ""
    traj = Traj3D()

    # Si l'on veut mesurer la distance entre le premier nucléotide et lui-même rajouté en bout de chaîne
    if appendend: 
        seq.append(traj[0])

    traj.compute(seq, rot_table)
    coords = traj.getTraj()
    score = fct_poids(coords)
    return score




