import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.tanh(x)

Nmin = 10 ** 2
Nmax = 10 ** 5
N_range = np.logspace(np.log10(Nmin), np.log10(Nmax), 1000, dtype=int)

a = 0
b = 1


random_numbers = [np.random.uniform(a, b, N) for N in N_range]
g_i = [func(rand) for rand in random_numbers] 

int_exact = np.log(np.cosh(b)) - np.log(np.cosh(a))
int_values = [(b - a) * np.mean(g) for g in g_i]
error = [np.abs(int_exact - int_val) for int_val in int_values]



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Left subplot
ax1.plot(N_range, int_values, label='int_values')
ax1.hlines(int_exact, Nmin, Nmax, label='int_exact', linestyle='--', color='red')
ax1.set_xlabel('N')
ax1.set_ylabel('Values')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()

# Right subplot
ax2.plot(N_range, error, label='error')
ax2.plot(N_range, [error[0] / np.sqrt(N/N_range[0]) for N in N_range], '--', label='error trend')
ax2.set_xlabel('N')
ax2.set_ylabel('Absolute Error')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.show()

