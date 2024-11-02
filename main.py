import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from monte_carlo_sim import MonteCarloSimulation
from parameters import *



def main():
    simulation = MonteCarloSimulation(N_photon, R, num_density, 
                               sigma, num_tracked_photons, photon_dist, **dist_params)
    simulation.simulate()

    mean_collisions = simulation.plot_coll_number_histogram()
    simulation.plot_energy_spectrum()
    simulation.plot_trajectories()


if __name__ == '__main__':
    main()