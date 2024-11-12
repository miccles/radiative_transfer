import warnings
import numpy as np
from monte_carlo_sim import MonteCarloSimulation
from parameters import *
from functions import TheoreticalDistributions
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# def sample_blackbody(theta_g):
#     def wien_peak_energy(theta_g):
#         temp_kelvin = theta_g * mec2_eV / kB
#         hv_peak = 2.4315 * 10 ** (-4) * temp_kelvin / mec2_eV  # returns peak energy in me * c^2 units
#         return hv_peak
#     while True:
#         print('______attempting _____')
#         print(f'wien_peak:{wien_peak_energy(theta_g)}')
#         en_rand = np.random.uniform(10 ** (-5), 10 ** 1) * wien_peak_energy(theta_g)
#         print(f'en_rand: {en_rand}')
#         max_loc = 2.82144 * theta_g
#         dist_max = TheoreticalDistributions(max_loc, 'blackbody', {'theta_g': theta_g}).probability_density()
#         print(f'dist_max: {dist_max}')
#         dist = TheoreticalDistributions(en_rand, 'blackbody', {'theta_g': theta_g}).probability_density()
#         print(f'dist: {dist}')
#         if np.random.random() < dist / dist_max:
#             return en_rand
        

# energy = sample_blackbody(0.01)
# print(energy)
# probability = TheoreticalDistributions(energy, 'blackbody', {'theta_g': 0.01}).probability_density()

import sys
import warnings

def enratio(x,mu):
    return 1 / (1 + x*(1 - mu))

def diff_cross_section(x, mu):
    return (3 / (16 * np.pi)) * enratio(x, mu) ** 2 * (enratio(x, mu) + 1 / enratio(x, mu) + mu ** 2 - 1)

def tot_cross_section(x):
    sigma_sigmaT = (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 
    return sigma_sigmaT


def prob_density(x, mu):
    return 2 * np.pi * diff_cross_section(x, mu) / tot_cross_section(x)


xvalues = [0.01, 0.1, 1, 5, 10, 100]

for x in xvalues:
    mu = np.linspace(-1, 1, 10000)
    thetas = np.arccos(mu)
    y = prob_density(x, mu)
    plt.plot(thetas, y, label=f'x = {x}')

plt.xlabel(r'$\theta_f$', fontsize=14)
plt.ylabel('Probability Density', fontsize=14)
plt.xscale('linear')
plt.yscale('log')
plt.legend()
plt.show()
