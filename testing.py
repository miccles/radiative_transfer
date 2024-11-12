import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from particles import generate_energy
from functions import TheoreticalDistributions

def en_final(en_ph_0, beta_el, en_el, theta_ph_fin, theta_el_in, theta_el_fin):
    return en_ph_0 * (1 - beta_el * np.cos(theta_el_in)) / (1 - beta_el * np.cos(theta_el_fin) + en_ph_0 / en_el * (1 - np.cos(theta_ph_fin)))

# Parameters for the distributions
photon_dist = 'blackbody'
photon_dist_params = {'theta_g': 0.01}

electron_dist = 'powerlaw'
electron_dist_params = {'alpha': 1.2, 'E_min': 0.15, 'E_max': 0.5}

# Generate blackbody distribution energies for en_ph_0
num_samples = 10000
photon_energies = [generate_energy(photon_dist, **photon_dist_params) for _ in range(num_samples)]

# Generate power law distribution energies for electrons
electron_energies = [generate_energy(electron_dist, **electron_dist_params) for _ in range(num_samples)]
gamma_electrons = [en + 1 for en in electron_energies]
beta_electrons = [np.sqrt(1 - 1 / gamma**2) for gamma in gamma_electrons]

# Generate random theta values
theta_ph_fin = np.arccos(2 * np.random.random(num_samples) - 1)
theta_el_in = np.arccos(2 * np.random.random(num_samples) - 1)
theta_el_fin = np.arccos(2 * np.random.random(num_samples) - 1)

# Calculate the final photon energies using en_final
final_photon_energies = [en_final(en_ph_0, beta_el, en_el, theta_ph_fin[i], theta_el_in[i], theta_el_fin[i])
                         for i, (en_ph_0, beta_el, en_el) in enumerate(zip(photon_energies, beta_electrons, electron_energies))]

# Plot histograms of the photon energies before and after the calculation
fig, ax = plt.subplots(figsize=(12, 6))

# Define logarithmic bins
log_bins = np.logspace(np.log10(min(photon_energies + final_photon_energies)), np.log10(max(photon_energies + final_photon_energies)), 40)

# Histogram of initial photon energies
hist, bin_edges = np.histogram(photon_energies, bins=log_bins)
bin_widths = np.diff(bin_edges)
normalized_hist = hist / (1 * bin_widths)
ax.step(bin_edges[:-1], normalized_hist, where='mid', alpha=0.6, color='b', label='Initial Photon Energies')

# Histogram of final photon energies
hist, bin_edges = np.histogram(final_photon_energies, bins=log_bins)
bin_widths = np.diff(bin_edges)
normalized_hist = hist / (1 * bin_widths)
ax.step(bin_edges[:-1], normalized_hist, where='mid', alpha=0.6, color='r', label='Final Photon Energies')

# Repeat the process N times
N = 10
for n in range(N):
    # Generate new random theta values
    theta_ph_fin = np.arccos(2 * np.random.random(num_samples) - 1)
    theta_el_in = np.arccos(2 * np.random.random(num_samples) - 1)
    theta_el_fin = np.arccos(2 * np.random.random(num_samples) - 1)

    # Calculate the final photon energies using en_final
    final_photon_energies = [en_final(en_ph_0, beta_el, en_el, theta_ph_fin[i], theta_el_in[i], theta_el_fin[i])
                             for i, (en_ph_0, beta_el, en_el) in enumerate(zip(final_photon_energies, beta_electrons, electron_energies))]

    # Add the histogram of the final photon energies
    hist, bin_edges = np.histogram(final_photon_energies, bins=log_bins)
    bin_widths = np.diff(bin_edges)
    normalized_hist = hist / (1 * bin_widths)
    ax.step(bin_edges[:-1], normalized_hist, where='mid', alpha=0.6, label=f'Final Photon Energies (Iteration {n+1})')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\frac{E_K}{m_e c^2}$')
ax.set_ylabel('Probability density')
ax.legend(loc='best')
ax.set_title('Photon Energy Distribution Before and After Calculation')

plt.tight_layout()
plt.show()

