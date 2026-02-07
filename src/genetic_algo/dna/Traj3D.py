import math

#For drawing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from genetic_algo.dna.RotTable import RotTable


class Traj3D:
    """Represents a 3D trajectory"""

    # Vertical translation (elevation) between two di-nucleotides
    __MATRIX_T = np.array(
        [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, -3.38/2],
        [0, 0, 0, 1]]
    )

    def __init__(self,want_to_plot=False):
        self.__Traj3D = []
        self.fig = plt.figure() if want_to_plot else None
        self.ax = plt.axes(projection='3d') if want_to_plot else None

    def getTraj(self) -> list:
        return self.__Traj3D

    def compute(self, dna_seq: str, rot_table: RotTable):
        N = len(dna_seq)
        self.__Traj3D = np.zeros((N,4))
        self.__Traj3D[0] = [0.0, 0.0, 0.0, 1.0]

        # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = np.eye(4)  # Identity matrix

        step_matrices = {}

        # On enregistre la position du premier nucléotide
        self.__Traj3D = [np.array([0.0, 0.0, 0.0, 1.0])]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                matrices_Rz[dinucleotide], matrices_Q[dinucleotide] = \
                    self.__compute_matrices(rot_table, dinucleotide)
                step_matrices[dinucleotide] = (
                    self.__MATRIX_T 
                    @ matrices_Rz 
                    @ matrices_Q 
                    @ matrices_Rz 
                    @ self.__MATRIX_T
                )

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix = total_matrix @ step_matrices[dinucleotide]
            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
            self.__Traj3D[i] = total_matrix[:,3]

    def __compute_matrices(self, rot_table: RotTable, dinucleotide: str):

        Omega = math.radians(rot_table.getTwist(dinucleotide))
        # Create rotation matrix of theta on Z axis
        matrices_Rz = \
            np.array([[math.cos(Omega/2), math.sin(Omega/2), 0, 0],
                        [-math.sin(Omega/2), math.cos(Omega/2), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])

        sigma = rot_table.getWedge(dinucleotide)
        delta = rot_table.getDirection(dinucleotide)
        alpha = math.radians(sigma)
        beta = math.radians(delta - 90)
        # Rotate of -beta on Z axis
        # Rotate of -alpha on X axis
        # Rotate of beta on Z axis
        matrices_Q = \
            np.array([[math.cos(-beta), math.sin(-beta), 0, 0],
                        [-math.sin(-beta), math.cos(-beta), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]]) \
            @ np.array([[1, 0, 0, 0],
                            [0, math.cos(-alpha), math.sin(-alpha), 0],
                            [0, -math.sin(-alpha), math.cos(-alpha), 0],
                            [0, 0, 0, 1]]) \
            @ np.array([[math.cos(beta), math.sin(beta), 0, 0],
                        [-math.sin(beta), math.cos(beta), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
        
        return matrices_Rz, matrices_Q

    def draw(self):
        xyz = np.array(self.__Traj3D)
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        self.ax.plot(x,y,z)
        self.ax.scatter(x[0],y[0],z[0],color="green")
        self.ax.scatter(x[-1],y[-1],z[-1],color="red")
        plt.show()

    def save_fig(self, filename: str):
        self.fig.savefig(filename)

    def save_coords(self, filename: str):
        with open(filename, 'w') as f:
            for i in range(len(self.__Traj3D)):
                f.write(f"{self.__Traj3D[i][0]},{self.__Traj3D[i][1]},{self.__Traj3D[i][2]}\n")
            f.close()