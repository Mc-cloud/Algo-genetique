import unittest
import random
import numpy as np
from src.genetic_algo.core.selection import * 
from unittest.mock import MagicMock

# --- The Test Class ---

class TestSelectionMethods(unittest.TestCase):

    def setUp(self):
        self.individuals = []
        for i in range(10):
            ind = MagicMock()
            ind.score = float(i)
            ind.__lt__ = lambda self, other: self.score < other.score
            self.individuals.append(ind)
        
        self.taux_selec = 0.5 
        self.expected_count = int(len(self.individuals) * self.taux_selec)
        self.expected_count_tournament = self.expected_count + int(0.1*len(self.individuals))

    def test_selection_elitiste(self):
        result = selection_elitiste(self.individuals, self.taux_selec)
        self.assertEqual(len(result), self.expected_count)
        scores = [ind.score for ind in result]
        self.assertEqual(scores, [0.0, 1.0, 2.0, 3.0, 4.0])

    def test_selection_tournament(self):
        result = selection_tournament(self.individuals, self.taux_selec)
        self.assertLessEqual(len(result), self.expected_count_tournament)

        for ind in result:
            self.assertIn(ind, self.individuals)

    def test_selection_roulette(self):
        result = selection_roulette(self.individuals, self.taux_selec)
        self.assertEqual(len(result), self.expected_count)

    def test_selection_roulette_exp(self):
        result = selection_roulette_exp(self.individuals, self.taux_selec)
        self.assertEqual(len(result), self.expected_count)

    def test_selection_rang_reel(self):
        result = selection_rang_reel(self.individuals, self.taux_selec)
        self.assertEqual(len(result), self.expected_count)

    def test_selection_rang_geometrique(self):
        result = selection_rang_geometrique(self.individuals, self.taux_selec)
        self.assertEqual(len(result), self.expected_count)

    def test_selection_wrapper(self):
        result_elitiste = selection(self.individuals, self.taux_selec, "elitiste")
        self.assertEqual(len(result_elitiste), self.expected_count)
        
        result_default = selection(self.individuals, self.taux_selec, "non_existent")
        self.assertEqual(len(result_default), self.expected_count)

if __name__ == '__main__':
    unittest.main()