import numpy as np

from parameters import *


def generate_energy(dist_type='normal', **kwargs):
    if dist_type == 'normal':
        mean = kwargs.get('mean')
        sigma = kwargs.get('sigma')
        return np.random.normal(mean, sigma)
    elif dist_type == 'uniform':
        min_val = kwargs.get('low')
        max_val = kwargs.get('high')
        return np.random.uniform(min_val, max_val)
    elif dist_type == 'powerlaw':
        alpha = kwargs.get('alpha')
        xmin = kwargs.get('xmin')
        xmax = kwargs.get('xmax')
        xi = np.random.random()
        if alpha == 1:
            return np.exp(np.log(xmax / xmin) * xi + np.log(xmin))
        else:
            return ((xmax ** (1 - alpha) - xmin ** (1 - alpha)) * xi + xmin ** (1 - alpha)) ** (1 / (1 - alpha))
    elif dist_type == 'blackbody':
        pass
    else:
        raise ValueError(f"Unsupported distribution type: {dist_type}")

class Particle:
    def __init__(self, mass, charge, energy_dist, **dist_params):
        self.mass = mass
        self.charge = charge
        self.energy = generate_energy(energy_dist, **dist_params) # energy in me * c^2 units

    
    def energy_keV(self):
        return self.energy * mec2_keV


class Photon(Particle):
    def __init__(self, energy_dist, track_trajectory=False, **dist_params):
        super().__init__(0, 0, energy_dist, **dist_params)
        self.energy = generate_energy(energy_dist, **dist_params)
        self.x = 0
        self.y = 0
        self.z = 0
        self.track_trajectory = track_trajectory
        self.trajectory = []
        self.collisions = 0

    def move(self, L, theta, phi):
        self.x += L * np.sin(theta) * np.cos(phi)
        self.y += L * np.sin(theta) * np.sin(phi)
        self.z += L * np.cos(theta)
        if self.track_trajectory:
            self.trajectory.append((self.x, self.y, self.z))
        self.collisions += 1
    

    def energy_to_wavelength(self):
        return lambda_db / self.energy


    def sigma_klein_nishina(self): # sigma / sigma_Thomson
        x = self.energy
        return (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 


    def compton_scatter(self, angle):  # returns the energy of the scattered photon in me * c^2 units # angle in radians
        return self.energy / (1 + self.energy * (1 - np.cos(angle)))
