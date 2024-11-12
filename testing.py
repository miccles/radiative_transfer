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

# 1. Open the output file in write mode
output_file = open('output.txt', 'w')

# 2. Redirect stdout (print statements) and stderr (errors) to the output file
sys.stdout = output_file
sys.stderr = output_file

# 3. Redirect warnings to the output file
warnings.simplefilter("always")  # Show all warnings
warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: \
    print(f"{filename}:{lineno}: {category.__name__}: {message}")

# Example prints, warnings, and errors
print("This is a print statement.")

warnings.warn("This is a warning.")

try:
    1 / 0  # This will raise an error
except ZeroDivisionError as e:
    print(f"Error: {e}")

# 4. Close the file after you're done
output_file.close()