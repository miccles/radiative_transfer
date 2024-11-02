import numpy as np
import matplotlib.pyplot as plt

from parameters import *

from particles import Photon

class MonteCarloSimulation:
    def __init__(
                self, N, R, n, sigma, num_tracked_photons, photon_dist,
                **dist_params
                ):
        self.N = N
        self.R = R
        self.n = n
        self.sigma = sigma
        self.num_tracked_photons = num_tracked_photons
        self.photons = [Photon(photon_dist, **dist_params) for _ in range(N)]
        self.tracked_photons = []
        self.select_random_photons(self.num_tracked_photons)



    def simulate(self):
        for photon in self.photons:
            while self.is_inside_sphere(photon):
                r1, r2, r3 = np.random.random(3)
                tau = -np.log(r1)
                L = self.calc_L_from_tau(tau)
                theta = np.arccos(2 * r2 - 1)
                phi = 2 * np.pi * r3
                photon.move(L, theta, phi)

    
    def calc_L_from_tau(self, tau):
        return tau / (self.n * self.sigma)

    def is_inside_sphere(self, photon):
        return photon.x**2 + photon.y**2 + photon.z**2 < self.R**2

    def get_collisions(self):
        return [photon.collisions for photon in self.photons]

    def select_random_photons(self, num=1):
        if num > self.N:
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
        energies = [photon.energy for photon in self.photons]
        plt.hist(energies, bins=20)
        plt.savefig('energy_spectrum.png')
        plt.close()