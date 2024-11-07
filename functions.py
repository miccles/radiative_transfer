import numpy as np
from scipy.special import kv


class TheoreticalDistributions: # Returns probability densities normalized to 1
    def __init__(self, dist_type, **kwargs):
        self.dist_type = dist_type
        self.kwargs = kwargs


    def probability_density(self, x):
        if self.dist_type == 'blackbody':
            return self.blackbody_distr(x, **self.kwargs)
        elif self.dist_type == 'maxwell_juttner':
            return self.maxwell_juttner_distr(x, **self.kwargs)
        elif self.dist_type == 'uniform':
            return self.uniform_distr(x, **self.kwargs)
        elif self.dist_type == 'powerlaw':
            return self.powerlaw_distr(x, **self.kwargs)
        elif self.dist_type == 'normal':
            return self.normal_distr(x, **self.kwargs)
        elif self.dist_type == 'monoenergetic':
            return self.monoenergetic_distr(x, **self.kwargs)
        else:
            raise ValueError(f"Unsupported distribution type: {self.dist_type}")

    def blackbody_distr(self, en, theta_g):
        return (15 / (np.pi ** 4 * theta_g ** 4)) * en ** 3 / (np.exp(en / theta_g) - 1)

    def maxwell_juttner_distr(self, gamma, theta):
        beta = np.sqrt(1 - 1 / gamma**2)
        return gamma ** 2 * beta * np.exp(-gamma / theta) / (theta * kv(2, 1 / theta))

    def uniform_distr(self, x, low, high):
        if low <= x <= high:
            return 1 / (high - low)
        else:
            return 0

    def powerlaw_distr(self, x, alpha, xmin, xmax):
        if xmin <= x <= xmax:
            return (alpha - 1) * (x / xmin) ** (-alpha) / (xmin ** (1 - alpha) - xmax ** (1 - alpha))
        else:
            return 0

    def normal_distr(self, x, mean, std):
        return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)

    def monoenergetic_distr(self, x, energy):
        std = 1e-6  # Very small standard deviation to simulate a delta function
        return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - energy) / std) ** 2)
