import numpy as np
import matplotlib.pyplot as plt
from parameters import *


from functions import TheoreticalDistributions


print(TheoreticalDistributions(1, 'blackbody', theta_g=1).probability_density())