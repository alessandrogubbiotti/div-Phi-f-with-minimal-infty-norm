import numpy as np
import matplotlib.pyplot as plt

def load_csv(name):
    return np.loadtxt(name, delimiter=",")

# Load data
u = load_csv("u.csv")
ux = load_csv("ux.csv")
uy = load_csv("uy.csv")
phix = load_csv("phi_x.csv")
phiy = load_csv("phi_y.csv")
infty_vals = load_csv("energies.csv")

N = u.shape[0]
X, Y = np.meshgrid(np.arange(N), np.arange(N))

# Helper function for edge plotting
def plot_edges(ax, vx, vy, title):
    norm = np.maximum(np.abs(vx), np.abs(vy))
    vmax = np.max(norm)

    for i in range(N):
        for j in range(N):
            c = np.abs(vx[i,j])
            ax.plot([i, (i+1)%N], [j, j], color=plt.cm.viridis(c / vmax), linewidth=2)

            c = np.abs(vy[i,j])
            ax.plot([i, i], [j, (j+1)%N], color=plt.cm.viridis(c / vmax), linewidth=2)

    ax.set_title(title)
    ax.set_aspect('equal')
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)

# ∇u + Φ
sumx = ux + phix
sumy = uy + phiy

# ∞-norm of original ∇u (for histogram baseline)
grad_u_infty = np.max(np.sqrt(ux**2 + uy**2))

# 1. Plot heatmap of u
plt.figure(figsize=(6, 5))
plt.imshow(u, origin='lower', cmap='inferno')
plt.colorbar(label='u(x)')
plt.title("Potential u(x)")
plt.tight_layout()
plt.savefig("plot_u.png")
plt.close()

# 2. ∇u edges
fig, ax = plt.subplots(figsize=(6, 6))
plot_edges(ax, ux, uy, "Gradient ∇u")
plt.tight_layout()
plt.savefig("plot_grad_u.png")
plt.close()

# 3. Φ edges
fig, ax = plt.subplots(figsize=(6, 6))
plot_edges(ax, phix, phiy, "Cycle field Φ")
plt.tight_layout()
plt.savefig("plot_phi.png")
plt.close()

# 4. ∇u + Φ edges
fig, ax = plt.subplots(figsize=(6, 6))
plot_edges(ax, sumx, sumy, "Improved field ∇u + Φ")
plt.tight_layout()
plt.savefig("plot_grad_u_plus_phi.png")
plt.close()

# 5. Histogram of ∞-norms
plt.figure(figsize=(7, 5))
plt.hist(infty_vals, bins=50, color='steelblue', alpha=0.7, label="||∇u + Φ||∞")
plt.axvline(grad_u_infty, color='red', linestyle='--', linewidth=2, label="||∇u||∞")
plt.xlabel("∞-norm")
plt.ylabel("Frequency")
plt.title("Histogram of ||∇u + Φ||∞")
plt.legend()
plt.tight_layout()
plt.savefig("hist_infty_norms.png")
plt.close()

