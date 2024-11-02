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


simulation.plot_coll_number_histogram()
simulation.plot_trajectories(num_photons=10)
