import numpy as np
import matplotlib.pyplot as plt
from parameters import *

mec2_eV = 511000  # Electron rest mass energy in eV

def sample_blackbody(theta_g):
    def wien_peak_energy(theta_g):
        temp_kelvin = theta_g * mec2_eV / kB
        en_peak = 2.4315 * 10 ** (-4) * temp_kelvin / mec2_eV  # returns peak energy in me * c^2 units
        return en_peak

    while True:
        en_rand = np.random.uniform(10 ** (-2), 10 ** 2) * wien_peak_energy(theta_g)
        dist = (15 / (np.pi ** 4 * theta_g ** 4)) * en_rand ** 3 / (np.exp(en_rand / theta_g) - 1)
        if np.random.random() < dist:
            return en_rand

def theoretical_blackbody(en, theta_g):
    return (15 / (np.pi ** 4 * theta_g ** 4)) * en ** 3 / (np.exp(en / theta_g) - 1)


# Parameters
theta_g = 1
num_samples = 10000

# Generate samples
samples = [sample_blackbody(theta_g) for _ in range(num_samples)]

# Plot histogram of samples
plt.hist(samples, bins=50, density=True, alpha=0.6, color='g', label='Sampled')

# Plot theoretical distribution
energy_values = np.linspace(0, max(samples), 1000)
theoretical_values = theoretical_blackbody(energy_values, theta_g)
plt.plot(energy_values, theoretical_values, 'r-', label='Theoretical')

plt.xlabel('Energy (in me * c^2 units)')
plt.ylabel('Probability density')
plt.legend()
plt.title('Blackbody Distribution')
plt.show()