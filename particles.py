import numpy as np
from scipy.optimize import newton
from parameters import *
from functions import TheoreticalDistributions 


def generate_energy(dist_type='normal', **kwargs): # RETURNS KINETIC ENERGY IN me * c^2 UNITS
    if dist_type == 'monoenergetic':
        return kwargs.get('energy')
    if dist_type == 'normal':
        mean = kwargs.get('mean')
        sigma = kwargs.get('sigma')
        return np.random.normal(mean, sigma)
    elif dist_type == 'uniform':
        min_val = kwargs.get('E_min')
        max_val = kwargs.get('E_max')
        return np.random.uniform(min_val, max_val)
    elif dist_type == 'powerlaw':
        alpha = kwargs.get('alpha')
        xmin = kwargs.get('E_min')
        xmax = kwargs.get('E_max')
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
    def func(gamma, theta):
        return gamma ** 3 - 2 * gamma ** 2 * theta - gamma + theta
    gamma_initial_guess = 1 + theta  # Initial guess for gamma
    gamma_max = newton(func, gamma_initial_guess, args=(theta,))
    while True:
        gamma_limit = 1 + 24.158 * theta ** 0.9478
        gamma_rand = np.random.uniform(1, gamma_limit)
        dist = TheoreticalDistributions.maxwell_juttner(gamma_rand - 1, theta)
        dist_max = TheoreticalDistributions.maxwell_juttner(gamma_max - 1, theta)
        if np.random.random() < dist / dist_max:
            return gamma_rand
        

def sample_blackbody(theta_g):
    def wien_peak_energy(theta_g):
        temp_kelvin = theta_g * mec2_eV / kB
        hv_peak = 2.4315 * 10 ** (-4) * temp_kelvin / mec2_eV  # returns peak energy in me * c^2 units
        return hv_peak
    while True:
        en_rand = np.random.uniform(10 ** (-5), 10 ** 1) * wien_peak_energy(theta_g)
        max_loc = 2.82144 * theta_g
        dist_max = TheoreticalDistributions(max_loc, 'blackbody', {'theta_g': theta_g}).probability_density()
        dist = TheoreticalDistributions(en_rand, 'blackbody', {'theta_g': theta_g}).probability_density()
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
    
    def lorentz_transform(self, electron, mu): # From lab frame to electron frame 
        gamma = 1 + electron.energy
        beta = np.sqrt(1 - 1 / gamma ** 2)
        en_photon = self.energy
        en_photon_prime = gamma * en_photon * (1 - beta * mu)
        return en_photon_prime

    def inverse_lorentz_transform(self, en_photon_prime, electron, mu): # From electron rest frame to lab frame
        gamma = 1 + electron.energy # Lab frame
        beta = np.sqrt(1 - 1 / gamma ** 2) # Lab frame
        en_photon = gamma * en_photon_prime * (1 + beta * mu) # Lab frame
        return en_photon


    def compton_energy_ratio(self, x, mu):
        return 1 / (1 + x * (1 - mu))

    def differential_cross_section(self, x, mu):
        enratio = self.compton_energy_ratio(x, mu)
        return (3 / (16 * np.pi)) * enratio ** 2 * (enratio + 1 / enratio + mu ** 2 - 1)

    def sigma_tot_klein_nishina(self, x): # sigma / sigma_Thomson
        sigma_sigmaT = (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 
        return sigma_sigmaT
    
    def angle_prob_density(self, x, mu):
        return 2 * np.pi * self.differential_cross_section(x, mu) / self.sigma_tot_klein_nishina(x)

    def sigma(self, x):
        return sigma_Thomson * self.sigma_tot_klein_nishina(x)
    

    def sample_angle(self, x): # generates scattering angle alpha from Klein-Nishina distribution
    # from alpha we get theta_f = theta_i - alpha
        while True:
            mu_rand = np.random.uniform(-1, 1)
            dist_max = self.angle_prob_density(x, 0)
            dist_value = self.angle_prob_density(x, mu_rand)
            if np.random.random() < dist_value / dist_max:
                theta_rand = np.arccos(mu_rand)
                return np.random.choice([theta_rand, -theta_rand])

