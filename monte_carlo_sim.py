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
        self.tracked_photons = []
        self.select_random_photons(self.num_tracked_photons)

    def simulate(self):
        for photon in self.photons:
            # Step 0: Initialize photons and move them to their initial positions
            r1, r2, r3 = np.random.random(3)
            tau = -np.log(r1)
            L = self.calc_L_from_tau(photon, tau)
            theta = np.arccos(2 * r2 - 1)
            phi = 2 * np.pi * r3
            photon.move(L, theta, phi)

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


    def calc_L_from_tau(self, photon, tau):
        sigma = photon.sigma()
        return tau / (self.n * sigma)
    

    def compton_scattering(self, photon):
        def en_final(en_ph_0, beta_el, en_el, theta_ph_fin, theta_el_in, theta_el_fin):
            return en_ph_0 * (1 - beta_el * np.cos(theta_el_in)) / (1 - beta_el * np.cos(theta_el_fin) + en_ph_0 / en_el * (1 - np.cos(theta_ph_fin)))
        
        # Generate an electron from the electron distribution
        electron = Particle(me, qe, self.electron_dist, **self.electron_dist_params)
        energy_electron = electron.energy
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
        plt.hist(collisions, bins=20)
        plt.vlines(np.mean(collisions), 0, 0.2 * N_photon, color='red', label='Mean')
        plt.vlines(max_tau ** 2, 0, 0.2 * N_photon, color='green', label=r'$\tau_{max}^2$')
        plt.legend(loc='best')
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
            # Plot the sampled photon energy histogram
            energies = [photon.energy for photon in self.photons]
            plt.hist(energies, bins=30, density=True, alpha=0.6, color='g', label='Sampled')

            # Plot the initial photon energy distribution
            initial_photon_energies = [generate_energy(self.photon_dist, **self.photon_dist_params) for _ in range(10000)]
            plt.hist(initial_photon_energies, bins=30, density=True, alpha=0.6, color='b', label='Initial Photon Distribution')

            # Plot the electron energy distribution
            electron_energies = [generate_energy(self.electron_dist, **self.electron_dist_params) for _ in range(10000)]
            plt.hist(electron_energies, bins=30, density=True, alpha=0.6, color='r', label='Electron Distribution')

            # Plot the photon initial distribution
            energy_values = np.linspace(0, max(energies), 1000)
            theoretical_values = TheoreticalDistributions(energy_values, self.photon_dist, self.photon_dist_params['theta_g'])
            plt.plot(energy_values, theoretical_values, 'r-', label='Theoretical')

            # Plot the electron initial distribution
            

            plt.xlabel(r'$\frac{E_K}{m_e c^2}$')
            plt.ylabel('Probability density')
            plt.legend(loc='best')
            plt.title('Energy Spectrum')