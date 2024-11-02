import numpy as np
import matplotlib.pyplot as plt

from particles import Photon

class MonteCarloSimulation:
    def __init__(
                self, N, R, n, sigma, photon_dist,
                **dist_params
                ):
        self.N = N
        self.R = R
        self.n = n
        self.sigma = sigma
        self.photons = [Photon(photon_dist, **dist_params) for _ in range(N)]

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

    def select_random_photon(self, num=1):
        if num > self.N:
            raise ValueError("Number of photons to select is greater than total number of photons")
        return np.random.choice(self.photons, num)

    def plot_coll_number_histogram(self):
        collisions = self.get_collisions()
        plt.hist(collisions, bins=20)
        plt.savefig('collisions_histogram.png')

    def plot_trajectories(self, num_photons):
            photons = self.select_random_photon(10)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            colors = plt.cm.jet(np.linspace(0, 1, len(photons)))
            
            for photon, color in zip(photons, colors):
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
            
            plt.savefig('trajectories.png')

    def plot_energy_spectrum(self):
        energies = [Photon.energy for photon in self.photons]
        plt.hist(energies, bins=20)
        plt.savefig('energy_spectrum.png')