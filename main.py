import numpy as np
import matplotlib.pyplot as plt

c = 1
me = 1
h = 1
lambda_db = h / (me * c) # de Broglie wavelength of electrons

mec2_keV = 510998.951 # electron rest mass energy in keV

N_photon = 1000

num_density = 1
sigma = 1
l_mean = 1 / (num_density * sigma)
R = 20 * l_mean

def calc_L_from_tau(tau, num_density, sigma):
    return tau / (num_density * sigma)


class Particle:
    def __init__(self, mass, charge, energy):
        self.mass = mass
        self.charge = charge
        self.energy = energy # energy in me * c^2 units

    
    def energy_keV(self):
        return self.energy * mec2_keV




class Photon(Particle):
    def __init__(self, energy):
        super().__init__(0, 0, energy)
        self.x = 0
        self.y = 0
        self.z = 0
        self.collisions = 0
        self.trajectory = [(self.x, self.y, self.z)]

    def move(self, L, theta, phi):
        self.x += L * np.sin(theta) * np.cos(phi)
        self.y += L * np.sin(theta) * np.sin(phi)
        self.z += L * np.cos(theta)
        self.collisions += 1
        self.trajectory.append((self.x, self.y, self.z))

    def energy_to_wavelength(self):
        return lambda_db / self.energy


    def sigma_klein_nishina(self): # sigma / sigma_Thomson
        x = self.energy
        return (3 / (8 * x)) * ((1 - 2 * (x + 1) / x ** 2) * np.log(1 + 2 * x) + 0.5 + 4 / x - 0.5 / (1 + 2 * x)**2) 


    def compton_scatter(self, angle):  # returns the energy of the scattered photon in me * c^2 units # angle in radians
        return self.energy / (1 + self.energy * (1 - np.cos(angle)))


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

    def select_random_photon(self):
        return np.random.choice(self.photons)

    def plot_trajectory(self, photon):
        trajectory = np.array(photon.trajectory)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Add a red dot at the origin
        ax.scatter(0, 0, 0, color='red', s=20)
        
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


random_photon = simulation.select_random_photon()
simulation.plot_trajectory(random_photon)