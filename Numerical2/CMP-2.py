import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 1.0                         # lattice constant
m_vals = np.arange(-3, 4)       # band index m
k_vals = np.linspace(-6*np.pi/a, 6*np.pi/a, 400)

'''
Q1: Periodicity = a
u_m = exp(i m 2pi x / a)
E_m(k) = (k + 2pi m / a)^2

'''
E_a = np.array([(k_vals + 2*np.pi*m/a)**2 for m in m_vals])

'''
Q2: Periodicity = 2a
G = pi m / a
E_m(k) = (k + pi m / a)^2
'''
E_2a = np.array([(k_vals + np.pi*m/a)**2 for m in m_vals])

# Plot Q1

plt.figure()
for i, m in enumerate(m_vals):
    plt.plot(k_vals, E_a[i], label=f"m = {m}")

plt.xlabel("k")
plt.ylabel("E")
plt.title("Free-electron bands (periodicity a)")
plt.legend(ncol=2)
plt.grid(True)
plt.show()


# Plot Q2

plt.figure()
for i, m in enumerate(m_vals):
    plt.plot(k_vals, E_2a[i], label=f"m = {m}")

plt.xlabel("k")
plt.ylabel("E")
plt.title("Free-electron bands (periodicity 2a)")
plt.legend(ncol=2)
plt.grid(True)
plt.show()
