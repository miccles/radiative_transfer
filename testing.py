import numpy as np
import matplotlib.pyplot as plt
from parameters import *




random_nums = np.random.random(1000)
tau = -np.log(random_nums)
mean_tau = np.mean(tau)
print(mean_tau)
plt.hist(tau, bins=20)
plt.vlines(mean_tau, 0, 0.2 * 1000, color='red', label='Mean')
plt.show()