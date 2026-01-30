from algo.fitness import *
import dna.RotTable as RotTable
import dna.Traj3D as Traj3D
import copy


## données
Table_rot = {
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
}
Nucleotids = [["AA","TT"],["AC","GT"],["AG", "CT"],["AT"],["CA","TG"],["CC","GG"],["CG"],["GA","TC"],["GC"],["TA"]]

base_seq = ''.join([line.rstrip('\n') for line in open("data/plasmid_8k.fasta")][1:]) 


##fonctions
def Rot_table_of(x) :
    D = copy.deepcopy(Table_rot)
    for i_couple in range(len(Nucleotids)):
        for nucleotid in Nucleotids[i_couple]:
            for k in range(2):
                D[nucleotid][k] = x[i_couple*2+k]
    return RotTable.RotTable(D)


def kernel(x, xp,sigma = 1,ell = 1):
    return sigma**2 * np.exp(-np.sum((x-xp)**2) / (2*ell**2))

##parametres
nb_space = 3
nb_cut,nbappend = 1,1
nb_ite = 20

## etat de depart
#Incertitudes = [( 0.06 ,  0.6),( 1.3  ,  5),(1.5  ,  3 ),( 1.1  ,  2),( 0.9  , 34 ),(  0.07 ,  2.1 ),(1.1  , 1.5),(0.9  ,  6),( 1.2  ,  1.275),[ 1.1  ,  2]]
Incertitudes = np.array([0.06,0.6,1.3,5,1.5,3,1.1,2,0.9,34,0.07,2.1,1.1,1.5,0.9,6,1.2,1.275,1.1,2])
#X = [np.array([(35.62 , 7.2),(34.4  , 1.1),(27.7  , 8.4),(31.5  , 2.6),(34.5  , 3.5),(33.67 , 2.1),(29.8  , 6.7),(36.9  , 5.3),(40    , 5 ),(36    , 0.9)])]
Space = []
X = [np.array([35.62,7.2,34.4,1.1,27.7,8.4,31.5,2.6,34.5,3.5,33.67,2.1,29.8,6.7,36.9,5.3,40,5,36,0.9])]
Y = [fitness(Rot_table_of(x),str_data,nbcuts=nb_cut,nbappend=nbappend) for x in X]

## opti bayesienne
for ite in range(nb_ite):
    K = np.array([[kernel(x, xp) for xp in X] for x in X])
    K += 1e-6 * np.eye(len(X))
    K_inv = np.linalg.inv(K)
    K_star = np.array([[kernel(x, xp) for xp in X] for x in Space])
    K_starstar = np.array([[kernel(x, xp) for xp in Space] for x in Space])
    mu_star = K_star @ K_inv @ Y
    cov_star = K_starstar - K_star @ K_inv @ K_star.T
    std_star = np.sqrt(np.diag(cov_star))