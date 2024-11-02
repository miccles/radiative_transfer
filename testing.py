import numpy as np
import matplotlib.pyplot as plt
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


electron_1 = Particle(me, -1, 1)
print(electron_1.energy)

photon_1 = Photon(2)
print(photon_1.energy)