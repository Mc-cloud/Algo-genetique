import numpy as np

class Individual:
    def __init__(self, base_rotation=None):
        if base_rotation is not None:
            self.chromosome = np.array(base_rotation, dtype = float)
        else :
            self.chromosome = np.zeros((16,3))
        self.fitness = 0.0


    def mutate(self, mutation_rate, sigma):
        mutation_mask = np.random.rand(16,3) < mutation_rate
        noise = np.random.normal(0, sigma, size = (16,3))
        self.chromosome += mutation_mask*noise

    