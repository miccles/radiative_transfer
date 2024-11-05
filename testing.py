import numpy as np
import matplotlib.pyplot as plt
from parameters import *


def blackbody_theor(en ,theta_g):
    return (15 / (np.pi ** 4 * theta_g ** 4)) * en ** 3 / (np.exp(en / theta_g) - 1)


energies = np.array([0.01, 0.1, 0.2, 0.32, 0.3])

energy_values = np.linspace(0, max(energies), 1000)
theoretical_values = blackbody_theor(energy_values, 0.01)


plt.plot(energy_values, theoretical_values, label='Theoretical blackbody distribution')
plt.show()