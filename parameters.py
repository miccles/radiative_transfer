c = 1
me = 1
h = 1

lambda_db = h / (me * c) # de Broglie wavelength of electrons

mec2_keV = 510998.951 # electron rest mass energy in keV

photon_dist = 'powerlaw'
alpha = 1.5
xmin = 0.01
xmax = 0.1
dist_params = {'alpha': alpha, 'xmin': xmin, 'xmax': xmax}


N_photon = 1000
num_tracked_photons = 0

num_density = 1
sigma = 1
l_mean_Thomson = 1 / (num_density * sigma)
max_tau = 20
R = max_tau * l_mean_Thomson