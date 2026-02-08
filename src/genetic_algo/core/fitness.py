import numpy as np
from genetic_algo.dna.RotTable import RotTable
from genetic_algo.dna.Traj3D import Traj3D


def dist_df(coords: list, nbappend = 1):
    """
    Takes as input a trajectory in the form of an array list 
    corresponding to the coordinates of the nodes, joining the first nbappend nodes
    at the end. The Euclidean distance between the
    initial nodes and the added nodes is then compared.“”"

    “”"
    Calculates the closure distance of a DNA trajectory.
    Measures the Euclidean distance between the first and last nodes
    of a trajectory to assess whether the DNA structure closes correctly.
    Args:
        coords: List of numpy arrays representing the 3D coordinates of the nodes of the trajectory
        nbappend: Number of nodes to compare between the start and end (default: 1)
    Returns:
        Total Euclidean distance between the nbappend first and last nodes
    Note:
        Used to check the circularity of DNA (circular plasmids)
    """
    start = coords[:nbappend]
    end = coords[-nbappend:]
    distsq = 0

    for i in range(nbappend):
        distsq+= np.linalg.norm(start[i]-end[i])**2

    return np.sqrt(distsq)

def dist_euclid(scores: list) -> float:
    """
    Calculates the Euclidean norm of a vector of scores.
    Args:
        scores: List of numerical values
    Returns:
        L2 (Euclidean) norm of the vector of scores
    Note:
        Used as a combination function to aggregate multiple scores
    """
    return np.linalg.norm(scores)


def fitness(rot_table: RotTable, seq: str, fct_poids = dist_df, nbappend = 2, nbcuts = 2, coup_combin = dist_euclid) -> float :
    """
    Fitness function to evaluate the quality of a DNA sequence.
    
    Calculates a closure score by testing the sequence at different cut-off points
    to evaluate its ability to form a stable circular structure.
    
    Args:
        rot_table: Rotation table for calculating the 3D trajectory of the DNA.
        seq: DNA sequence to be evaluated (character string).
        fct_weight: Function for calculating the weight/score for a trajectory (default: dist_df).
        nbappend: Number of nodes to add at the end to test closure (default: 2)
        nbcuts: Number of cut points to test (default: 2)
        coup_combin: Function to combine the scores of the different cuts (default: dist_euclid)
    
    Returns:
        Fitness score (the lower the score, the better the closure)
    
    Note:
        The method tests several circular rotations of the sequence to evaluate
        the overall stability of the circular DNA structure.
    """
    nbases = len(seq)
    assert nbases >= nbappend
    traj = Traj3D()

    def eval_une_coupure(seq: str, nbappend: int, indcut: int):
        """
        Evaluates the closure score for a given cut site.
        Args:
            seq: DNA sequence
            nbappend: Number of nodes to add
            indcut: Cut site index
        
        Returns:
            Closure score for this cut
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