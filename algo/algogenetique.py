import dna.RotTable as RotTable
import dna.Traj3D as Traj3D
import numpy as np
import random
import fitness
import select
from json import load as json_load

str_data = 'AACTGTCAGCTACCGATCATCTAGCTCTATATCGCGCATTAGCAGC'
Rot_data = {
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

class Individu:
    def __init__(self, Table_rot):
        self.Rot_table = RotTable(Table_rot)
        self.score = self.fit()


    def __add__(self, other): #fonction permettant d'acoupler deux individus
        global str_data
        s_score,o_score = self.score,other.score
        Table ={}
        for XY in self.Rot_table:
            Ls = self.Rot_table[XY]
            Lo = other.Rot_table[XY]
            for i in range(3):
                alpha = 0.5  # implementer alpha en fct des scores
                Ls[i] = alpha*Ls[i] +(1-alpha)*Lo[i]
            Table[XY] = Ls
        return Individu(Table)
    
    def mutation(self,mutrate,sigma):
        global Rot_data
        if 0<=mutrate <=1:
            Table = self.Rot_table
            for XY in Table:
                L = Table[XY]
                for i in range(3):
                    if random.random() < mutrate :
                        tamp = L[i] + random.normalvariate(0,sigma*L[i+3])
                        L[i] = max(Rot_data[XY][i]-Rot_data[XY][i+3],min(tamp,Rot_data[XY][i]+Rot_data[XY][i+3]))
                Table[XY] = L
            self.Rot_table = Table

    def fit(self) -> float: #Renvoie le score de l'individu
        return fitness.fitness(self.Rot_table,str_data)

    def __lt__(self,other):
        return self.score<other.score
    
class Population :
    def __init__(self,filename : str,dna_seq: str, nb_individus, nb_generations,taux_selec,selection_type : str):
        self.rot_table = json_load(open(filename))
        self.Individus = self.New_pop()
        self.dna = dna_seq
        self.taux_selec = taux_selec
        self.nb_individus = nb_individus
        self.nb_generations = nb_generations
        self.best_table = self.Start_Algo_gene()
    
    def New_pop(self):
        def New_individu():
            Table_rot = self.rot_table
            for XY in Table_rot:
                L = Table_rot[XY]
                for i in range(3):
                        L[i] = random.uniform(-Rot_data[XY][i+3],Rot_data[XY][i+3])
                Table_rot[XY] = L
            return Individu(Table_rot)
            
        L = [New_individu() for _ in range (self.nb_individus)]
        self.Individus = L


    def Start_Algo_gene(self):
        1
