import numpy as np

from parameters import *

class Particle:
    def __init__(self, mass, charge, energy):
        self.mass = mass
        self.charge = charge
        self.energy = energy # energy in me * c^2 units

    
    def energy_keV(self):
        return self.energy * mec2_keV




class Photon(Particle):
    def __init__(self, energy):
        super().__init__(0, 0, energy)
        self.x = 0
        self.y = 0
        self.z = 0
        self.collisions = 0
        self.trajectory = [(self.x, self.y, self.z)]

    def move(self, L, theta, phi):
        self.x += L * np.sin(theta) * np.cos(phi)
        self.y += L * np.sin(theta) * np.sin(phi)
        self.z += L * np.cos(theta)
        self.collisions += 1
        self.trajectory.append((self.x, self.y, self.z))

    def energy_to_wavelength(self):
        return lambda_db / self.energy


    def sigma_klein_nishina(self): # sigma / sigma_Thomson
        x = self.energy
        return (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 


    def compton_scatter(self, angle):  # returns the energy of the scattered photon in me * c^2 units # angle in radians
        return self.energy / (1 + self.energy * (1 - np.cos(angle)))
