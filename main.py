import numpy as np
import matplotlib.pyplot as plt

from particles import Particle, Photon
from parameters import *



N_photon = 1000

num_density = 1
sigma = 1
l_mean = 1 / (num_density * sigma)
R = 20 * l_mean

def calc_L_from_tau(tau, num_density, sigma):
    return tau / (num_density * sigma)




class MonteCarloSimulation:
    def __init__(self, N, R, n, sigma, energy):
        self.N = N
        self.R = R
        self.n = n
        self.sigma = sigma
        self.photons = [Photon(energy) for _ in range(N)]

    def simulate(self):
        for photon in self.photons:
            while self.is_inside_sphere(photon):
                r1, r2, r3 = np.random.random(3)
                tau = -np.log(r1)
                L = tau / (self.n * self.sigma)
                theta = np.arccos(2 * r2 - 1)
                phi = 2 * np.pi * r3
                photon.move(L, theta, phi)

    def is_inside_sphere(self, photon):
        return photon.x**2 + photon.y**2 + photon.z**2 < self.R**2

    def get_collisions(self):
        return [photon.collisions for photon in self.photons]

    def select_random_photon(self, num=1):
        if num > self.N:
            raise ValueError("Number of photons to select is greater than total number of photons")
        return np.random.choice(self.photons, num)

    def plot_trajectory(self, photons):
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
            
            plt.show()


simulation = MonteCarloSimulation(N_photon, R, num_density, sigma, 0.01)
simulation.simulate()
collisions = simulation.get_collisions()

plt.hist(collisions, bins=20)
plt.show()


random_photons = simulation.select_random_photon(10)
simulation.plot_trajectory(random_photons)