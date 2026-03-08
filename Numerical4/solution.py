import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# Construct Hamiltonian and compute eigenvalues (band energy)
# -----------------------------------------------------------

def compute_bands(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk):

    kmin = -bz_index*np.pi/lattice_a
    kmax =  bz_index*np.pi/lattice_a

    kgrid = np.linspace(kmin, kmax, Nk)

    energies = np.zeros((nbasis, Nk))

    G0 = 2*np.pi/lattice_a
    Glist = np.arange(-nbasis//2, nbasis//2 + 1)[:nbasis] * G0

    for ik, k in enumerate(kgrid):

        Hmat = np.zeros((nbasis, nbasis))

        for i in range(nbasis):
            for j in range(nbasis):

                if i == j:
                    Hmat[i,i] = 0.5*(k - Glist[i])**2
                else:
                    order = abs(i-j)
                    Hmat[i,j] = Vcoeff.get(order,0)

        eigvals = np.linalg.eigvalsh(Hmat)
        energies[:,ik] = eigvals

    return kgrid, energies


# -----------------------------------------------------------
# Plot band structure
# -----------------------------------------------------------

def draw_band_structure(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk):

    k, bands = compute_bands(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk)

    plt.figure(figsize=(7,6))

    for n in range(bands.shape[0]):
        plt.plot(k, bands[n], lw=1)

    # Brillouin zone boundaries
    if bz_index == 1:
        edges = [-np.pi/lattice_a, np.pi/lattice_a]
    else:
        edges = [-2*np.pi/lattice_a, -np.pi/lattice_a,
                 np.pi/lattice_a,  2*np.pi/lattice_a]

    for e in edges:
        plt.axvline(e, ls="--", color="gray", lw=0.7)

    plt.xlabel("k")
    plt.ylabel("Energy (eV)")
    plt.title(f"Band structure ({bz_index} BZ)")
    plt.show()


# -----------------------------------------------------------
# Density of States using Gaussian broadening
# -----------------------------------------------------------

def compute_DOS(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk):

    k, bands = compute_bands(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk)

    Emin = bands.min()
    Emax = bands.max()

    Egrid = np.linspace(Emin, Emax, 800)

    sigma = 0.1
    DOS = np.zeros_like(Egrid)

    for band in bands:
        for Ek in band:
            DOS += np.exp(-(Egrid-Ek)**2/(2*sigma**2))

    return Egrid, DOS


def draw_DOS(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk):

    E, D = compute_DOS(A0, lattice_a, nbasis, bz_index, Vcoeff, Nk)

    plt.figure()
    plt.plot(E, D)
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    plt.title("Density of States")
    plt.show()


# -----------------------------------------------------------
# Fermi Energy (1 electron per unit cell)
# -----------------------------------------------------------

def fermi_energy(bands):

    first_band = bands[0]
    Ef = 0.5*(first_band.max()+first_band.min())

    return Ef


# -----------------------------------------------------------
# Parameters
# -----------------------------------------------------------

A = 3
a = 1
Nk = 2000

Vcoeff = {1:A/2, 2:A/4, 3:A/6}

# -----------------------------------------------------------
# 1. Band structure
# -----------------------------------------------------------

draw_band_structure(A, a, 5, 1, Vcoeff, Nk)
draw_band_structure(A, a, 5, 2, Vcoeff, Nk)

# -----------------------------------------------------------
# 2. Density of states
# -----------------------------------------------------------

draw_DOS(A, a, 5, 2, Vcoeff, Nk)


# -----------------------------------------------------------
# 3. Convergence with Fourier components
# -----------------------------------------------------------

basis_sets = [
    {1:A/2},
    {1:A/2,2:A/4},
    {1:A/2,2:A/4,3:A/6}
]

for V in basis_sets:
    draw_band_structure(A,a,5,2,V,Nk)


# -----------------------------------------------------------
# 4. Fermi energy
# -----------------------------------------------------------

k,bands = compute_bands(A,a,5,1,Vcoeff,Nk)

Ef = fermi_energy(bands)
print("Fermi energy =",Ef)

plt.figure()
for b in bands:
    plt.plot(k,b)

plt.axhline(Ef,color="red",ls="--")
plt.title("Fermi level")
plt.show()


# convergence with k-mesh
Nk_list = [100,200,400,800,1500,3000]
Ef_list = []

for n in Nk_list:
    k,b = compute_bands(A,a,5,1,Vcoeff,n)
    Ef_list.append(fermi_energy(b))

plt.figure()
plt.plot(Nk_list,Ef_list,'o-')
plt.xlabel("Number of k points")
plt.ylabel("Fermi energy")
plt.show()


# -----------------------------------------------------------
# 5. Modified potential
# -----------------------------------------------------------

Vnew = {1:A,2:A/2,3:A/4,4:A/6}

draw_band_structure(A,2*a,5,2,Vnew,Nk)
draw_DOS(A,2*a,5,2,Vnew,Nk)


# -----------------------------------------------------------
# 6. Total energy per unit cell
# -----------------------------------------------------------

Avals = np.linspace(0,5,12)

Ecell = []

for Aval in Avals:

    Vtemp = {1:Aval/2,2:Aval/4,3:Aval/6}

    k,b = compute_bands(Aval,a,5,1,Vtemp,Nk)

    Ef = fermi_energy(b)

    occ = b[b<=Ef]

    Etot = occ.sum()/len(k)

    Ecell.append(Etot)

plt.figure()
plt.plot(Avals,Ecell,'o-')
plt.xlabel("A (eV)")
plt.ylabel("Energy per unit cell")
plt.title("Total energy vs potential amplitude")
plt.show()