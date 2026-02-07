# tests/test_algogenetique.py

import unittest
import src.genetic_algo.core.algogenetique  as alg


str_data = 'AACTGTCAGCTACCGATCATCTAGCTCTATATCGCGCATTAGCAGCCAGCATCGACATCGTAGCTCACGCGATATCCGATCGTAGCGCTGCGAGCGCTGCTAGCTAGCTAGTCGATGCATGCTAGCTACGATGCAT'
Rot_data_place = "src/genetic_algo/dna/table.json"

class TestAlgoGenetique(unittest.TestCase):
    def setUp(self):
        _ = alg.AlgoGenetique(Rot_data_place,str_data,2,1,0.5,'')
        Rot_data = alg.Rot_data
        self.Individu1 = alg.Individu(Rot_data)
        self.Individu2 = alg.Individu(Rot_data)

    def test_mutation(self):
        Rot_data = alg.Rot_data
        New_individu = self.Individu1
        New_individu.mutation(1,10)
        D = New_individu.Rot_table.rot_table
        for XY in D :
            for i in range(3):
                self.assertAlmostEqual(D[XY][i], Rot_data[XY][i], delta=Rot_data[XY][i+3]+1e-6)
    
    def test_AlgoGenetique(self):
        Best,_ ,_= alg.AlgoGenetique(Rot_data_place,str_data,8,4,0.5,'')
        self.assertGreaterEqual(self.Individu1.score,Best[-1].score)

if __name__ == "__main__":
    unittest.main()

