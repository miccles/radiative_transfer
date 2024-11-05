import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from monte_carlo_sim import MonteCarloSimulation
from parameters import *



def main():
    photon_dist_params = get_particle_params('photon', photon_dist)
    electron_dist_params = get_particle_params('electron', electron_dist)
    simulation = MonteCarloSimulation(
        N_photon, num_tracked_photons,
        R, num_density, 
        photon_dist=photon_dist, photon_dist_params=photon_dist_params,
        electron_dist=electron_dist, electron_dist_params=electron_dist_params)
    simulation.simulate()

    mean_collisions = simulation.plot_coll_number_histogram()
    simulation.plot_energy_spectrum()
    simulation.plot_trajectories()


if __name__ == '__main__':
    main()


