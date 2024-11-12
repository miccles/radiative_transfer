import numpy as np
import matplotlib.pyplot as plt

from parameters import *
from functions import TheoreticalDistributions
from particles import Photon, Particle, generate_energy


class MonteCarloSimulation:
    def __init__(
                self, num_photons, num_tracked_photons,
                R, n, 
                photon_dist, photon_dist_params, 
                electron_dist, electron_dist_params
                ):
        self.num_photons = num_photons
        self.num_tracked_photons = num_tracked_photons
        self.R = R
        self.n = n
        self.num_tracked_photons = num_tracked_photons
        self.photon_dist = photon_dist
        self.photon_dist_params = photon_dist_params
        self.electron_dist = electron_dist
        self.electron_dist_params = electron_dist_params
        self.photons = [Photon(self.photon_dist, **self.photon_dist_params) for _ in range(num_photons)]
        # for photon in self.photons:
        #     if photon.energy < 10 ** (-10):
        #         print(f'Warning, energy = {photon.energy}')
        #     else:
        #         print(f'Photon Initilization, energy = {photon.energy}')
        self.tracked_photons = []
        self.select_random_photons(self.num_tracked_photons)
        self.cross_sections = []

    def simulate(self):
        tracker = 0
        for photon in self.photons:
            # tracker += 1
            # print('---------------------------------')
            # print(f'Simulating photon {tracker}...')
            #Step 0: Initialize photons and move them to their initial positions
            r1, r2, r3 = np.random.random(3)
            tau = -np.log(r1)
            L = self.calc_L_from_tau(photon, tau)
            theta = np.arccos(2 * r2 - 1)
            phi = 2 * np.pi * r3
            photon.move(L, theta, phi)
            self.cross_sections.append(photon.sigma())

            # Apply Compton scattering and update photon energy
            while self.is_inside_sphere(photon):
                en_photon_f, theta_photon_f = self.compton_scattering(photon)
                photon.energy = en_photon_f

                # Generate two random numbers for L and phi
                r1, r2 = np.random.random(2)
                tau = -np.log(r1)
                L = self.calc_L_from_tau(photon, tau)
                phi = 2 * np.pi * r2

                # Propagate the photon
                photon.move(L, theta_photon_f, phi)
            #print(f'Photon {tracker} escaped. Ncollisions = {photon.collisions}, Energy = {photon.energy}')

    def calc_L_from_tau(self, photon, tau):
        sigma = photon.sigma()
        return tau / (self.n * sigma)
    

    def compton_scattering(self, photon):
        def en_final(en_ph_0, beta_el, en_el, theta_ph_fin, theta_el_in, theta_el_fin):
            # print('..............')
            # print(f'en_ph_0: {en_ph_0}')
            # print(f'beta_el: {beta_el}')
            # print(f'en_el: {en_el}')
            # print(f'theta_ph_fin: {theta_ph_fin}')
            # print(f'theta_el_in: {theta_el_in}')
            # print(f'theta_el_fin: {theta_el_fin}')
            return en_ph_0 * (1 - beta_el * np.cos(theta_el_in)) / (1 - beta_el * np.cos(theta_el_fin) + en_ph_0 / en_el * (1 - np.cos(theta_ph_fin)))
        
        # Generate an electron from the electron distribution
        electron = Particle(me, qe, self.electron_dist, **self.electron_dist_params)
        energy_electron = electron.energy
        # if energy_electron < 10 ** (-10):
        #     print(f'Electron initilization, potential overflow')
        # if photon.energy < 10 ** (-100):
        #     print(f'WARNING: extremely low photon energy')
        # elif photon.energy < 10 ** (-10):
        #     print(f'Warning: photon energy={photon.energy}')
        gamma_electron = energy_electron + 1
        beta_electron = np.sqrt(1 - 1 / gamma_electron**2)

        # Generate angles from an isotropic distribution
        r1, r2, r3 = np.random.random(3)
        theta_el_in = np.arccos(2 * r1 - 1)
        theta_el_fin = np.arccos(2 * r2 - 1)
        theta_photon_f = np.arccos(2 * r3 - 1)

        # Calculate the final energy of the photon
        en_photon_f = en_final(photon.energy, beta_electron, energy_electron, theta_photon_f, theta_el_in, theta_el_fin)
        return en_photon_f, theta_photon_f

    def is_inside_sphere(self, photon):
        return photon.x**2 + photon.y**2 + photon.z**2 < self.R**2

    def get_collisions(self):
        return [photon.collisions for photon in self.photons]

    def select_random_photons(self, num=1):
        if num > self.num_photons:
            raise ValueError("Number of photons to select is greater than total number of photons")
        selected_photons = np.random.choice(self.photons, num, replace=False)
        for photon in selected_photons:
            photon.track_trajectory = True
            photon.trajectory = [(photon.x, photon.y, photon.z)]
            self.tracked_photons.append(photon)

    def plot_coll_number_histogram(self):
        collisions = self.get_collisions()
        cross_sections = self.cross_sections 
        Ntot = len(collisions)

        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot the collision number histogram with logarithmic bins, normalized to 1
        log_bins = np.logspace(np.log10(min(collisions)), np.log10(max(collisions)), 30)
        hist, bin_edges = np.histogram(collisions, bins=log_bins)
        bin_widths = np.diff(bin_edges)
        normalized_hist = hist / (Ntot * bin_widths)
        axs[0].bar(bin_edges[:-1], normalized_hist, width=bin_widths, align='edge', alpha=0.6, color='r')
        axs[0].vlines(np.mean(collisions), 0, max(normalized_hist), color='red', label='Mean')
        axs[0].vlines(max_tau ** 2, 0, max(normalized_hist), color='green', label=r'$\tau_{max}^2$')
        axs[0].set_xscale('log')
        axs[0].set_yscale('log')
        axs[0].set_xlabel('Number of Collisions')
        axs[0].set_ylabel('Probability density')
        axs[0].legend(loc='best')
        axs[0].set_title('Collision Number Histogram')

        # Plot the cross-section histogram with linear bins, normalized to 1
        hist, bin_edges = np.histogram(cross_sections, bins=40, density=True)
        axs[1].hist(cross_sections, bins=40, density=True, alpha=0.6, color='orange')
        axs[1].set_xlabel(r'$\sigma/\sigma_T$')
        axs[1].set_ylabel('Probability density')
        axs[1].set_title('Cross Section Histogram')

        plt.tight_layout()
        plt.savefig('collisions_histogram.png')
        plt.close()

        return np.mean(collisions)

    def plot_trajectories(self):
        if not self.tracked_photons:
            return
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        colors = plt.cm.jet(np.linspace(0, 1, len(self.tracked_photons)))

        for photon, color in zip(self.tracked_photons, colors):
            if photon.trajectory is not None:
                trajectory = np.array(photon.trajectory)
                ax.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2], color=color)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Add a red dot at the origin
        ax.scatter(0, 0, 0, color='red', s=100)

        # Plot a transparent sphere of radius R
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = self.R * np.cos(u) * np.sin(v)
        y = self.R * np.sin(u) * np.sin(v)
        z = self.R * np.cos(v)
        ax.plot_wireframe(x, y, z, color="gray", alpha=0.3)

        plt.show()

    def plot_energy_spectrum(self):
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        # Plot the initial photon energy distribution + theoretical distribution
        initial_photon_energies = [generate_energy(self.photon_dist, **self.photon_dist_params) for _ in range(1000)]
        log_bins = np.logspace(np.log10(min(initial_photon_energies)), np.log10(max(initial_photon_energies)), 30)
        hist, bin_edges = np.histogram(initial_photon_energies, bins=log_bins)
        bin_widths = np.diff(bin_edges)
        normalized_hist = hist / (1000 * bin_widths)
        axs[0, 0].bar(bin_edges[:-1], normalized_hist, width=bin_widths, align='edge', alpha=0.6, color='b', label=r'$f_{\gamma}(E_i)$')
        energy_values = np.linspace(min(initial_photon_energies), max(initial_photon_energies), 1000)
        theoretical_values = TheoreticalDistributions(energy_values, self.photon_dist, self.photon_dist_params).probability_density()
        axs[0, 0].plot(energy_values, theoretical_values, 'r-', label=r'$f_{\gamma, \text{th}}(E_i)$')
        axs[0, 0].set_xlabel(r'$\frac{E_K}{m_e c^2}$')
        axs[0, 0].set_ylabel('Probability density')
        axs[0, 0].set_xscale('log')
        axs[0, 0].set_yscale('log')
        #axs[0, 0].set_xlim(0.9 * min(initial_photon_energies), 1.1 * max(initial_photon_energies))
        axs[0, 0].legend(loc='best')
        axs[0, 0].set_title('Initial Photon Energy Distribution')
        
        

        # Plot the electron energy distribution + theoretical distribution
        electron_energies = [generate_energy(self.electron_dist, **self.electron_dist_params) for _ in range(1000)]
        log_bins = np.logspace(np.log10(min(electron_energies)), np.log10(max(electron_energies)), 30)
        hist, bin_edges = np.histogram(electron_energies, bins=log_bins)
        bin_widths = np.diff(bin_edges)
        normalized_hist = hist / (1000 * bin_widths)
        axs[0, 1].bar(bin_edges[:-1], normalized_hist, width=bin_widths, align='edge', alpha=0.6, color='b', label=r'$f_{e}(E)$')
        energy_values = np.linspace(min(electron_energies), max(electron_energies), 1000)
        theoretical_values = TheoreticalDistributions(energy_values, self.electron_dist, self.electron_dist_params).probability_density()
        axs[0, 1].plot(energy_values, theoretical_values, 'r-', label=r'$f_{e, \text{th}}(E)$')
        axs[0, 1].set_xlabel(r'$\frac{E_K}{m_e c^2}$')
        axs[0, 1].set_ylabel('Probability density')
        axs[0, 1].set_xscale('log')
        axs[0, 1].set_yscale('log')
        #axs[0, 1].set_xlim(0.9 * min(electron_energies), 1.1 * max(electron_energies))
        axs[0, 1].legend(loc='best')
        axs[0, 1].set_title('Electron Energy Distribution')

        # Plot the sampled photon energy histogram
        energies = [photon.energy for photon in self.photons]
        energy_threshold = 10 ** (-5)
        energies = [energy for energy in energies if energy > energy_threshold]
        log_bins = np.logspace(np.log10(min(energies)), np.log10(max(energies)), 30)
        hist, bin_edges = np.histogram(energies, bins=log_bins)
        bin_widths = np.diff(bin_edges)
        normalized_hist = hist / (self.num_photons * bin_widths)
        axs[1, 0].bar(bin_edges[:-1], normalized_hist, width=bin_widths, align='edge', alpha=0.6, color='g', label=r'$f_{\gamma}(E_f)$')
        axs[1, 0].set_xlabel(r'$\frac{E_K}{m_e c^2}$')
        axs[1, 0].set_ylabel('Probability density')
        axs[1, 0].set_xscale('log')
        axs[1, 0].set_yscale('log')
        axs[1, 0].legend(loc='best')
        axs[1, 0].set_title('Final Photon Energy Spectrum')

        # Hide the empty subplot (bottom right)
        axs[1, 1].axis('off')

        plt.tight_layout()
        plt.savefig('energy_spectrum.png')
        plt.close()