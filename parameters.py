c = 1
me = 1
qe = -1
h = 1

lambda_db = h / (me * c) # de Broglie wavelength of electrons

mec2_eV = 510998.951 # electron rest mass energy in eV
kB = 8.617333262145 * 10 ** (-5) # Boltzmann constant in eV / K

N_photon = 5
num_tracked_photons = 0

num_density = 1
sigma_Thomson = 1
l_mean_Thomson = 1 / (num_density * sigma_Thomson)
max_tau = 20
R = max_tau * l_mean_Thomson

photon_dist = 'blackbody'
electron_dist = 'powerlaw'

# Photon distribution parameters
photon_dist_params_dict = {
    'monoenergetic': {'energy': 0.01},
    'normal': {'mean': 0.01, 'std': 0.1},
    'uniform': {'E_min': 0.01, 'E_max': 10},
    'powerlaw': {'alpha': 1, 'E_min': 0.01, 'E_max': 10},
    'blackbody': {'theta_g': 0.01},
    'maxwell_juttner': {'theta': 1, 'gamma_max': 10}
}

# Electron distribution parameters
electron_dist_params_dict = {
    'monoenergetic': {'energy': 1.1},
    'normal': {'mean': 1.1, 'std': 0.1},
    'uniform': {'E_min': 1, 'E_max': 1.2},
    'powerlaw': {'alpha': 1.2, 'E_min': 0.15, 'E_max': 0.5},
    'blackbody': {'theta_g': 1},
    'maxwell_juttner': {'theta': 1, 'gamma_max': 10}
}

def get_particle_params(particle_type, dist_type):
    if particle_type == 'photon':
        return photon_dist_params_dict.get(dist_type, {})
    elif particle_type == 'electron':
        return electron_dist_params_dict.get(dist_type, {})


