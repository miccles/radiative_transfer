import numpy as np
from scipy.special import kv



def blackbody_distr(en ,theta_g):
    return (15 / (np.pi ** 4 * theta_g ** 4)) * en ** 3 / (np.exp(en / theta_g) - 1)


def maxwell_juttner_distr(gamma, theta):
    beta = np.sqrt(1 - 1 / gamma**2)
    return gamma ** 2 * beta * np.exp(-gamma / theta) /(theta * kv(2, 1 / theta)) 