import warnings
import sys

# Open the output file in write mode
output_file = open('output.txt', 'w')
# Redirect stdout (print statements) and stderr (errors) to the output file
sys.stdout = output_file
sys.stderr = output_file

# Redirect warnings to the output file
warnings.simplefilter("always")  # Show all warnings
warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: \
    print(f"{filename}:{lineno}: {category.__name__}: {message}")

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


# Close the file after the script finishes
output_file.close()