import numpy as np
import matplotlib.pyplot as plt

hbar=1.0
m=1.0
a=1.0
G=2*np.pi/a

U={1:1.2,2:0.6,3:0.3}

def solve_nfe(k,U,n_basis=10):
    n=np.arange(-n_basis,n_basis+1)
    size=len(n)
    H=np.zeros((size,size),dtype=complex)
    for i in range(size):
        for j in range(size):
            if i==j:
                H[i,j]=(hbar**2*(k-n[i]*G)**2)/(2*m)
            else:
                diff=abs(n[i]-n[j])
                H[i,j]=U.get(diff,0)
    return np.linalg.eigvalsh(H)

k_vec=np.linspace(-3.5*np.pi/a,3.5*np.pi/a,2000)
dk=k_vec[1]-k_vec[0]

energies=[]
for k in k_vec:
    energies.append(solve_nfe(k,U))
energies=np.array(energies)

fig,ax=plt.subplots(figsize=(12,7))
colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728']

for i in range(4):
    ax.plot(k_vec/(np.pi/a),energies[:,i],
            label=f'Band {i+1}',
            color=colors[i],lw=2.5)

bz=[-3,-2,-1,0,1,2,3]
for b in bz:
    if b!=0:
        ax.axvline(b,color='red',ls='--',alpha=0.3)

plt.title(r"$\mathbf{Extended\ Zone\ Scheme:\ 1st,\ 2nd,\ and\ 3rd\ BZs}$")
plt.xlabel(r"Wavevector $k$ [units of $\pi/a$]")
plt.ylabel(r"Energy $E$")
plt.xticks(bz,
           [r"$-3\pi/a$",r"$-2\pi/a$",r"$-\pi/a$",r"$0$",
            r"$\pi/a$",r"$2\pi/a$",r"$3\pi/a$"])
plt.legend(loc='upper left')
plt.grid(True,linestyle=':',alpha=0.6)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12,6))

for i in range(3):
    slope=np.gradient(energies[:,i],dk)
    plt.plot(k_vec/(np.pi/a),slope,
             label=f'Slope Band {i+1}',
             color=colors[i],lw=2)

for b in bz:
    plt.plot(b,0,'ko',markersize=8)
    plt.axvline(b,color='red',ls='--',alpha=0.3)

plt.axhline(0,color='black',lw=1.2)
plt.title(r"$\mathbf{Zero\ Slope\ at\ All\ BZ\ Boundaries}$")
plt.xlabel(r"Wavevector $k$ [units of $\pi/a$]")
plt.ylabel(r"Slope $dE/dk$")
plt.xticks(bz,
           [r"$-3\pi/a$",r"$-2\pi/a$",r"$-\pi/a$",r"$0$",
            r"$\pi/a$",r"$2\pi/a$",r"$3\pi/a$"])
plt.ylim(-15,15)
plt.legend()
plt.grid(True,alpha=0.3)
plt.tight_layout()
plt.show()