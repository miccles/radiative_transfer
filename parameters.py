c = 1
me = 1
h = 1

lambda_db = h / (me * c) # de Broglie wavelength of electrons

mec2_keV = 510998.951 # electron rest mass energy in keV

photon_dist = 'monoenergetic'
photon_energy = 0.01
photon_dist_params = {'energy': photon_energy}

electorn_dist = 'monoenergetic'
electron_energy = 1.1
electorn_dist_params = {'energy': electron_energy}


N_photon = 1000
num_tracked_photons = 0

num_density = 1
sigma_Thomson = 1
l_mean_Thomson = 1 / (num_density * sigma_Thomson)
max_tau = 20
R = max_tau * l_mean_Thomson