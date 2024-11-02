import numpy as np
import matplotlib.pyplot as plt

from monte_carlo_sim import MonteCarloSimulation
from parameters import *


simulation = MonteCarloSimulation(N_photon, R, num_density, sigma, photon_dist, alpha=alpha, xmin=xmin, xmax=xmax
                                  )
simulation.simulate()


simulation.plot_coll_number_histogram()
simulation.plot_trajectories(num_photons=10)
simulation.plot_energy_spectrum()