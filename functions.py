import numpy as np
from scipy.special import kv


class TheoreticalDistributions: # Returns probability densities normalized to 1
    def __init__(self, energy, dist_type, **kwargs):
        self.energy = energy     # kinetic energy in me * c^2 units - equialent to gamma - 1
        self.dist_type = dist_type
        self.kwargs = kwargs


    def probability_density(self):
        if self.dist_type == 'blackbody':
            return self.blackbody()
        elif self.dist_type == 'maxwell_juttner':
            return self.maxwell_juttner()
        elif self.dist_type == 'uniform':
            return self.uniform()
        elif self.dist_type == 'powerlaw':
            return self.powerlaw()
        elif self.dist_type == 'normal':
            return self.normal()
        elif self.dist_type == 'monoenergetic':
            return self.monoenergetic()
        else:
            raise ValueError(f"Unsupported distribution type: {self.dist_type}")

    def blackbody(self):
        en = self.energy
        theta_g = self.kwargs.get('theta_g')
        return (15 / (np.pi ** 4 * theta_g ** 4)) * en ** 3 / (np.exp(en / theta_g) - 1)

    def maxwell_juttner(self):
        gamma = self.energy + 1
        theta = self.kwargs.get('theta')
        beta = np.sqrt(1 - 1 / gamma**2)
        return gamma ** 2 * beta * np.exp(-gamma / theta) / (theta * kv(2, 1 / theta))

    def uniform(self):
        low = self.kwargs.get('E_min')
        high = self.kwargs.get('E_max')
        return 1 / (high - low)

    def powerlaw(self):
        x = self.energy
        alpha = self.kwargs.get('alpha')
        xmin = self.kwargs.get('E_min')
        xmax = self.kwargs.get('E_max')
        if xmin <= x <= xmax:
            return (alpha - 1) * (x / xmin) ** (-alpha) / (xmin ** (1 - alpha) - xmax ** (1 - alpha))
        else:
            return 0

    def normal(self):
        x = self.energy
        mean = self.kwargs.get('mean')
        std = self.kwargs.get('std')
        return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)

    def monoenergetic(self):
        pass
    #     x = self.energy
    #     std = 1e-6  # Very small standard deviation to simulate a delta function
    #     return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - energy) / std) ** 2)
