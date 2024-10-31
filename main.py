import numpy as np
import matplotlib.pyplot as plt

from monte_carlo_sim import MonteCarloSimulation
from parameters import *



N_photon = 1000

num_density = 1
sigma = 1
l_mean = 1 / (num_density * sigma)
R = 20 * l_mean


simulation = MonteCarloSimulation(N_photon, R, num_density, sigma, 0.01)
simulation.simulate()
collisions = simulation.get_collisions()

plt.hist(collisions, bins=20)
plt.show()


random_photons = simulation.select_random_photon(10)
simulation.plot_trajectory(random_photons)