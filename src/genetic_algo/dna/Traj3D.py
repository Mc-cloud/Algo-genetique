import math

#For drawing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from genetic_algo.dna.RotTable import RotTable

from numba import jit

@jit(nopython=True)
def fast_compute_loop(encoded_seq, matrices_db):
    """
    Compiled machine code version of the trajectory loop.
    """
    N = len(encoded_seq)
    # Initialize trajectory array (N, 4)
    # We use flattened index for the last dimension to keep it simple in Numba
    traj = np.zeros((N, 4))
    
    # Starting point: [0, 0, 0, 1]
    traj[0, 0] = 0.0
    traj[0, 1] = 0.0
    traj[0, 2] = 0.0
    traj[0, 3] = 1.0

    # Initialize cumulative matrix as Identity (4x4)
    total_matrix = np.eye(4)

    # Loop over the sequence
    for i in range(1, N):
        # Calculate dinucleotide index: previous_base * 4 + current_base
        # (A=0, C=1, G=2, T=3) -> 'AC' = 0*4 + 1 = 1
        idx = encoded_seq[i-1] * 4 + encoded_seq[i]
        
        # Retrieve the pre-computed matrix for this dinucleotide
        step_matrix = matrices_db[idx]
        
        # Matrix multiplication: Total = Total @ Step
        # Numba optimizes this dot product heavily
        total_matrix = total_matrix @ step_matrix
        
        # Store the position (last column of the transformation matrix)
        traj[i, 0] = total_matrix[0, 3]
        traj[i, 1] = total_matrix[1, 3]
        traj[i, 2] = total_matrix[2, 3]
        traj[i, 3] = total_matrix[3, 3]

    return traj


class Traj3D:
    """Represents a 3D trajectory"""

    __NUCL_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    __DINUC_KEYS = [n1+n2 for n1 in "ACGT" for n2 in "ACGT"]

    # Vertical translation (elevation) between two di-nucleotides
    __MATRIX_T = np.array(
        [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, -3.38/2],
        [0, 0, 0, 1]],
        dtype = float
    )

    def __init__(self,want_to_plot=False):
        self.__Traj3D = []
        self.fig = plt.figure() if want_to_plot else None
        self.ax = plt.axes(projection='3d') if want_to_plot else None

    def getTraj(self) -> list:
        return self.__Traj3D

    def compute(self, dna_seq: str, rot_table: RotTable):
        """
        Computes the trajectory. Uses Numba if available for high performance.
        """
        # 1. Convert DNA string to integer array (A=0, C=1, G=2, T=3)
        # This is much faster for the computer to process than strings
        # We use a simple list comprehension or numpy map
        try:
            encoded_seq = np.array([self.__NUCL_MAP[s] for s in dna_seq], dtype=np.int32)
        except KeyError:
            # Fallback for unexpected characters like 'N' (treat as 'A' or handle error)
            encoded_seq = np.zeros(len(dna_seq), dtype=np.int32) 

        # 2. Pre-compute the 16 possible step matrices
        # We do this in Python because trigonometry is fast enough for just 16 items
        matrices_db = np.zeros((16, 4, 4))
        
        for i, dinuc in enumerate(self.__DINUC_KEYS):
            Rz, Q = self.__compute_matrices(rot_table, dinuc)
            matrices_db[i] = (
                self.__MATRIX_T 
                @ Rz 
                @ Q 
                @ Rz 
                @ self.__MATRIX_T
            )

        # 3. Run the high-performance loop
        self.__Traj3D = fast_compute_loop(encoded_seq, matrices_db)

    def __compute_matrices(self, rot_table: RotTable, dinucleotide: str):


        Omega = math.radians(rot_table.getTwist(dinucleotide))
        # Create rotation matrix of theta on Z axis
        

        sigma = rot_table.getWedge(dinucleotide)
        delta = rot_table.getDirection(dinucleotide)
        alpha = math.radians(sigma)
        beta = math.radians(delta - 90)

        cos_omega = math.cos(Omega/2)
        sin_omega = math.sin(Omega/2)
        cos_beta = math.cos(beta)
        sin_beta = math.sin(beta)
        cos_nbeta = math.cos(-beta)
        sin_nbeta = math.sin(-beta)
        cos_nalpha = math.cos(-alpha)
        sin_nalpha = math.sin(-alpha)

        matrices_Rz = np.array([
            [cos_omega, sin_omega, 0, 0],
            [-sin_omega, cos_omega, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ], dtype=float)

        # Constructing Q manually (or via multiplication)
        Q_rot_z1 = np.array([
            [cos_nbeta, sin_nbeta, 0, 0],
            [-sin_nbeta, cos_nbeta, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ], dtype=float)

        Q_rot_x = np.array([
            [1, 0, 0, 0],
            [0, cos_nalpha, sin_nalpha, 0],
            [0, -sin_nalpha, cos_nalpha, 0],
            [0, 0, 0, 1]
        ], dtype=float)

        Q_rot_z2 = np.array([
            [cos_beta, sin_beta, 0, 0],
            [-sin_beta, cos_beta, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ], dtype=float)

        matrices_Q = Q_rot_z1 @ Q_rot_x @ Q_rot_z2
        
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