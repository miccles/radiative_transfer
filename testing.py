import numpy as np
import matplotlib.pyplot as plt


N = 2000  # Number of random numbers

random_numbers_x = np.random.uniform(0, 1, N)
random_numbers_y = 1 - random_numbers_x

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Left subplot: Histogram of random_numbers
ax1.hist(random_numbers_x, bins=20)
ax1.set_title('Histogram of cosines')
ax1.set_xlabel('Value')
ax1.set_ylabel('Frequency')

# Right subplot: Histogram of cosines
ax2.hist(random_numbers_y, bins=20)
ax2.set_title('Histogram of phi angles')
ax2.set_xlabel('Value')
ax2.set_ylabel('Frequency')

plt.show()