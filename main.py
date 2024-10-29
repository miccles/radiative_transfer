import numpy as np
import matplotlib.pyplot as plt

c = 1
me = 1
h = 1
lambda_db = h / (me * c) # de Broglie wavelength of electrons

mec2_keV = 510998.951 # electron rest mass energy in keV


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


    def energy_to_wavelength(self):
        return lambda_db / self.energy


    def sigma_klein_nishina(self): # sigma / sigma_Thomson
        x = self.energy
        return (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 


    def compton_scatter(self, angle):  # returns the energy of the scattered photon in me * c^2 units # angle in radians
        return self.energy / (1 + self.energy * (1 - np.cos(angle)))



xrange = np.logspace(-3, 1, 1000)
