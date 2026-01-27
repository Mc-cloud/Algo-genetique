import unittest
import numpy as np
import os
# Importing real classes from your package structure
import dna.RotTable as RotTable
import dna.Traj3D as Traj3D
from algo.fitness import dist_df, dist_euclid, fitness, fitness_basic

class TestFitnessReal(unittest.TestCase):

    def setUp(self):
        """Set up real data for testing."""
        # Standard table for testing
        self.rot_table = RotTable({
            "AA": [35.62, 7.2, -154, 0.06, 0.6, 0],
            "AC": [34.4, 1.1, 143, 1.3, 5, 0],
            "AG": [27.7, 8.4, 2, 1.5, 3, 0],
            "AT": [31.5, 2.6, 0, 1.1, 2, 0],
            "CA": [34.5, 3.5, -64, 0.9, 34, 0],
            "CC": [33.67, 2.1, -57, 0.07, 2.1, 0],
            "CG": [29.8, 6.7, 0, 1.1, 1.5, 0],
            "CT": [27.7, 8.4, -2, 1.5, 3, 0],
            "GA": [36.9, 5.3, 120, 0.9, 6, 0],
            "GC": [40, 5, 180, 1.2, 1.275, 0],
            "GG": [33.67, 2.1, 57, 0.07, 2.1, 0],
            "GT": [34.4, 1.1, -143, 1.3, 5, 0],
            "TA": [36, 0.9, 0, 1.1, 2, 0],
            "TC": [36.9, 5.3, -120, 0.9, 6, 0],
            "TG": [34.5, 3.5, 64, 0.9, 34, 0],
            "TT": [35.62, 7.2, 154, 0.06, 0.6, 0]
        })
        self.test_seq = "AACTGTCAGCTACCGATCATCTAGCTCTATATCGCGCATTAGCAGC"

    def test_dist_df(self):
        """Verifies dist_df calculation with real numpy arrays."""
        coords = [
            np.array([0, 0, 0]), 
            np.array([1, 0, 0]), 
            np.array([2, 0, 0]), 
            np.array([0, 1, 0])
        ]
        # nbappend=1 compares [0,0,0] and [0,1,0] -> dist is 1.0
        self.assertEqual(dist_df(coords, nbappend=1), 1.0)

    def test_fitness_type(self):
        score = fitness(self.rot_table, self.test_seq, nbappend=2, nbcuts=1)
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0)

    def test_fitness_consistency(self):
        res1 = fitness_basic(self.rot_table, self.test_seq)
        res2 = fitness(self.rot_table, self.test_seq, nbcuts=0)
        self.assertEqual(res1, res2)

    def test_sequence_length_error(self):
        with self.assertRaises(AssertionError):
            fitness(self.rot_table, "A", nbappend=10)

if __name__ == '__main__':
    unittest.main()