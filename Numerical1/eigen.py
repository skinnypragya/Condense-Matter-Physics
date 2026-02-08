import numpy as np
import matplotlib.pyplot as plt
n = 100
I = 0
mat = np.zeros((n, n))

k = 1
for i in range(1, k+1):
    e = -5
    b = 1
    m = n
    k = 0
    while mat[k][k] != 0:
        k+=1
    for j in range(k, k+m):
        mat[j][j] = e
        if j+1 <= k+m-1:
            mat[j][j+1] = -b
        if j <= k+m-1 and j-1 >= k:
            mat[j][j-1] = -b
    for k in range(n//2-1):
        mat[k][k+4] = I
        mat[k+4][k] = I

# print(mat)

eigenvalues, eigenvectors = np.linalg.eig(mat)

eigenvalues = np.sort(eigenvalues)

for i in range(len(eigenvalues)//2):
    plt.hlines(eigenvalues[i], 0, 2, colors='b')
    plt.hlines(eigenvalues[-(i+1)], 0, 2, colors='r', alpha=0.5)
plt.title(f'Eigen Values for N = {n} atoms')
plt.xlabel('Index')
plt.ylabel('Eigen Value')
plt.text(1, eigenvalues[-1]+0.1, 'Red: Unoccupied Levels\nBlue: Occupied Fermi Levels', fontsize=10, verticalalignment='top', horizontalalignment='right')
plt.show()

# plot frequency of eigen values under specific bins
bw = float(input("Enter bin width for histogram: "))
bins = np.arange(np.min(eigenvalues), np.max(eigenvalues)+bw, bw)
plt.hist(eigenvalues, bins, color='blue')
plt.title(f'Histogram of Eigen Values with bin width {bw}')
plt.xlabel('Eigen Value')
plt.ylabel('Frequency')
plt.show()



def fermi_level(eigenvalues):
    N = len(eigenvalues)
    # each level holds 2 electrons → half filled
    return eigenvalues[N//2 - 1], eigenvalues[N//2]

def plot_eigenstates(matrix):
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)

    idx_min = 0
    idx_max = -1
    idx_homo = len(eigenvalues)//2 - 1
    idx_lumo = len(eigenvalues)//2

    states = {
        "Minimum Energy": idx_min,
        "HOMO": idx_homo,
        "LUMO": idx_lumo,
        "Maximum Energy": idx_max
    }

    x = np.arange(1, len(eigenvalues)+1)

    plt.figure(figsize=(10,6))
    for label, idx in states.items():
        plt.plot(x, eigenvectors[:, idx],
                 marker='o',
                 label=f"{label} (E={eigenvalues[idx]:.2f} eV)")

    plt.xlabel("Atomic site index")
    plt.ylabel("Eigenstate amplitude")
    plt.title(f"Eigenstates for 1D Chain (N={len(eigenvalues)})")
    plt.legend()
    plt.grid(True)
    plt.show()

plot_eigenstates(mat)