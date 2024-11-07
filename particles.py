import numpy as np
from scipy.optimize import newton
from parameters import *
from functions import * 


def generate_energy(dist_type='normal', **kwargs): # RETURNS KINETIC ENERGY IN me * c^2 UNITS
    if dist_type == 'monoenergetic':
        return kwargs.get('energy')
    if dist_type == 'normal':
        mean = kwargs.get('mean')
        sigma = kwargs.get('sigma')
        return np.random.normal(mean, sigma)
    elif dist_type == 'uniform':
        min_val = kwargs.get('Emin')
        max_val = kwargs.get('Emax')
        return np.random.uniform(min_val, max_val)
    elif dist_type == 'powerlaw':
        alpha = kwargs.get('alpha')
        xmin = kwargs.get('Emin')
        xmax = kwargs.get('Emax')
        xi = np.random.random()
        if alpha == 1:
            return np.exp(np.log(xmax / xmin) * xi + np.log(xmin))
        else:
            return ((xmax ** (1 - alpha) - xmin ** (1 - alpha)) * xi + xmin ** (1 - alpha)) ** (1 / (1 - alpha))
    elif dist_type == 'blackbody':
        theta_g = kwargs.get('theta_g')
        return sample_blackbody(theta_g)
    elif dist_type == 'maxwell_juttner':
        theta = kwargs.get('theta')
        gamma_max = kwargs.get('gamma_max')
        gamma_particle = sample_maxwell_juttner(theta, gamma_max)
        return gamma_particle
    else:
        raise ValueError(f"Unsupported distribution type: {dist_type}")


def sample_maxwell_juttner(theta, gamma_max):
    gamma_initial_guess = 1 + theta  # Initial guess for gamma
    gamma_max = newton(maxwell_juttner_distr, gamma_initial_guess, args=(theta,))
    while True:
        gamma_rand = np.random.uniform(1, gamma_max)
        dist = maxwell_juttner_distr(gamma_rand, theta)
        dist_max = maxwell_juttner_distr(gamma_max, theta)
        if np.random.random() < dist / dist_max:
            return gamma_rand
        

def sample_blackbody(theta_g):
    def wien_peak_energy(theta_g):
        temp_kelvin = theta_g * mec2_eV / kB
        hv_peak = 2.4315 * 10 ** (-4) * temp_kelvin / mec2_eV  # returns peak energy in me * c^2 units
        return hv_peak
    while True:
        en_rand = np.random.uniform(10 ** (-5), 10 ** 3) * wien_peak_energy(theta_g)
        max_loc = 2.82144 * theta_g
        dist_max = blackbody_distr(max_loc, theta_g)
        dist = blackbody_distr(en_rand, theta_g)
        if np.random.random() < dist / dist_max:
            return en_rand

class Particle:
    def __init__(self, mass, charge, energy_dist, **dist_params):
        self.mass = mass
        self.charge = charge
        self.energy = generate_energy(energy_dist, **dist_params) # energy in me * c^2 units

    
    def energy_eV(self):
        return self.energy * mec2_eV # energy in eV units


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


    def sigma(self):
        return sigma_Thomson * self.sigma_klein_nishina()

    def compton_scatter(self, angle):  # returns the energy of the scattered photon in me * c^2 units # angle in radians
        return self.energy / (1 + self.energy * (1 - np.cos(angle)))
