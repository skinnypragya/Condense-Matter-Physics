import numpy as np
import matplotlib.pyplot as plt


def ladder_model():

    print("="*60)
    print("Two Leg Ladder Band Structure")
    print("="*60)

    a = 1
    k = np.linspace(-np.pi/a, np.pi/a, 500)

    t = -1
    ratios = [0.5, 1, 2, 3, 4, 5]

    fig, ax = plt.subplots(2, 3, figsize=(14,9))
    ax = ax.ravel()

    for i, r in enumerate(ratios):

        tp = r*abs(t)

        e1 = 2*t*np.cos(k*a) + tp
        e2 = 2*t*np.cos(k*a) - tp

        ax[i].plot(k*a/np.pi, e1, lw=2)
        ax[i].plot(k*a/np.pi, e2, lw=2)

        ax[i].axhline(0, ls="--", alpha=0.5)

        gap = np.min(e1) - np.max(e2)

        phase = "INSULATOR" if gap>0 else "METAL"
        color = "red" if gap>0 else "green"

        ax[i].set_title(f"t'/t = {r}  |  Gap = {gap:.2f}  |  {phase}", color=color)
        ax[i].set_xlabel("ka/pi")
        ax[i].set_ylabel("Energy")
        ax[i].grid(alpha=0.3)

    plt.suptitle("Two Leg Ladder Dispersion")
    plt.tight_layout()
    plt.savefig("ladder_bands.png", dpi=150)
    plt.show()


def diagonal_hopping():

    print("\n"+"="*60)
    print("Diagonal Hopping Study")
    print("="*60)

    a = 1
    k = np.linspace(-np.pi/a, np.pi/a, 500)

    t = -1
    tp = 2

    tpp_list = [0, 0.3, 0.5, 0.7, 0.9]

    fig, ax = plt.subplots(2, len(tpp_list), figsize=(18,9))

    for i, tpp in enumerate(tpp_list):

        base = 2*t*np.cos(k*a)

        r = tp + tpp*np.cos(k*a)
        im = tpp*np.sin(k*a)

        mag = np.sqrt(r**2 + im**2)

        ep = base + mag
        em = base - mag

        gap = np.min(ep) - np.max(em)
        phase = "INSULATOR" if gap>0 else "METAL"

        ax[0,i].plot(k*a/np.pi, ep)
        ax[0,i].plot(k*a/np.pi, em)
        ax[0,i].set_title(f"Parallel  t''={tpp}  {phase}")
        ax[0,i].grid()

        off = np.abs(tp + 2*tpp*np.cos(k*a))

        ep = base + off
        em = base - off

        gap = np.min(ep) - np.max(em)
        phase = "INSULATOR" if gap>0 else "METAL"

        ax[1,i].plot(k*a/np.pi, ep)
        ax[1,i].plot(k*a/np.pi, em)
        ax[1,i].set_title(f"Crossed  t''={tpp}  {phase}")
        ax[1,i].grid()

    plt.suptitle("Effect of Diagonal Hopping")
    plt.tight_layout()
    plt.savefig("diagonal_effect.png", dpi=150)
    plt.show()


def square_lattice_nn():

    print("\n"+"="*60)
    print("2D Square Lattice (Nearest Neighbor)")
    print("="*60)

    a = 1
    t = -2

    N = 200

    kx = np.linspace(-np.pi/a, np.pi/a, N)
    ky = np.linspace(-np.pi/a, np.pi/a, N)

    KX, KY = np.meshgrid(kx, ky)

    energy = 2*t*(np.cos(KX*a) + np.cos(KY*a))

    fig, ax = plt.subplots(figsize=(6,6))

    c = ax.contourf(KX, KY, energy, 50, cmap="RdBu_r")
    plt.colorbar(c)

    ax.set_title("Energy Contour (NN)")

    plt.tight_layout()
    plt.savefig("square_nn.png")
    plt.show()


def square_lattice_nnn():

    print("\n"+"="*60)
    print("Square Lattice with NNN hopping")
    print("="*60)

    a = 1
    t = -2

    N = 200

    kx = np.linspace(-np.pi/a, np.pi/a, N)
    ky = np.linspace(-np.pi/a, np.pi/a, N)

    KX, KY = np.meshgrid(kx, ky)

    for tp in [-1, 1]:

        energy = 2*t*(np.cos(KX*a)+np.cos(KY*a)) + 4*tp*(np.cos(KX*a)*np.cos(KY*a))

        fig, ax = plt.subplots(figsize=(6,6))

        c = ax.contourf(KX, KY, energy, 50, cmap="RdBu_r")
        plt.colorbar(c)

        ax.set_title(f"NNN hopping  t' = {tp}")

        plt.tight_layout()
        plt.savefig(f"square_nnn_{tp}.png")
        plt.show()


if __name__ == "__main__":

    print("="*60)
    print("Running All Calculations")
    print("="*60)

    ladder_model()
    diagonal_hopping()
    square_lattice_nn()
    square_lattice_nnn()

    print("\nFinished Successfully")