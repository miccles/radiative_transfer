import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.tanh(x)

Nmin = 10 ** 2
Nmax = 10 ** 8
N_range = np.logspace(np.log10(Nmin), np.log10(Nmax), 1000, dtype=int)

a = 0
b = 1


random_numbers = [np.random.uniform(a, b, N) for N in N_range]
g_i = [func(rand) for rand in random_numbers] 

int_exact = np.log(np.cosh(b)) - np.log(np.cosh(a))
int_values = [(b - a) * np.mean(g) for g in g_i]
error = [np.abs(int_exact - int_val) for int_val in int_values]


plt.plot(N_range, error)
plt.plot(N_range, [error[0] / np.sqrt(N/N_range[0]) for N in N_range], '--')
plt.xlabel('N')
plt.ylabel('Absolute Error')
plt.xscale('log')
plt.yscale('log')
plt.show()
