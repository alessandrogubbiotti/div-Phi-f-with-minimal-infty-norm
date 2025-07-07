import numpy as np
import matplotlib.pyplot as plt

def load_field(filename):
    return np.loadtxt(filename, delimiter=',')

def compute_grad_u(u):
    N = u.shape[0]
    grad_u_x = np.zeros((N, N))
    grad_u_y = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            grad_u_x[i, j] = u[i, (j + 1) % N] - u[i, j]
            grad_u_y[i, j] = u[(i + 1) % N, j] - u[i, j]
    return grad_u_x, grad_u_y

def compute_vertex_vectors(x_edges, y_edges):
    N = x_edges.shape[0]
    h = 1.0 / N

    u_vert = np.zeros((N, N))
    v_vert = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            # Outgoing edges from vertex (i,j):
            # Horizontal edge from (i,j) to (i, j+1)
            x_out = x_edges[i, j]
            # Vertical edge from (i,j) to (i+1, j)
            y_out = y_edges[i, j]

            # Average over outgoing edges divided by lattice spacing
            # since grad ~ delta u / h
            u_vert[i, j] = x_out / h
            v_vert[i, j] = y_out / h

    return u_vert, v_vert

def plot_vertex_quiver(u, u_x, u_y, phi_x, phi_y, title):
    grad_u_x, grad_u_y = compute_grad_u(u)
    sum_x = grad_u_x + phi_x
    sum_y = grad_u_y + phi_y

    # Compute vertex vectors (mean of outgoing edges divided by h)
    u_vert_grad = compute_vertex_vectors(grad_u_x, grad_u_y)
    u_vert_phi = compute_vertex_vectors(phi_x, phi_y)
    u_vert_sum = compute_vertex_vectors(sum_x, sum_y)

    N = u.shape[0]
    X, Y = np.meshgrid(np.arange(N), np.arange(N))

    fig, axs = plt.subplots(1, 3, figsize=(18,6))

    axs[0].quiver(X, Y, u_vert_grad[0], u_vert_grad[1], color='r')
    axs[0].set_title(r'Vertex vectors of $\nabla u$')
    axs[0].set_aspect('equal')
    axs[0].invert_yaxis()
    axs[0].set_xticks([])
    axs[0].set_yticks([])

    axs[1].quiver(X, Y, u_vert_phi[0], u_vert_phi[1], color='g')
    axs[1].set_title(r'Vertex vectors of $\Phi$')
    axs[1].set_aspect('equal')
    axs[1].invert_yaxis()
    axs[1].set_xticks([])
    axs[1].set_yticks([])

    axs[2].quiver(X, Y, u_vert_sum[0], u_vert_sum[1], color='b')
    axs[2].set_title(r'Vertex vectors of $\nabla u + \Phi$')
    axs[2].set_aspect('equal')
    axs[2].invert_yaxis()
    axs[2].set_xticks([])
    axs[2].set_yticks([])

    plt.suptitle(title)
    plt.show()

def main():
    u = load_field("u.csv")
    phi_x = load_field("phi_x.csv")
    phi_y = load_field("phi_y.csv")

    plot_vertex_quiver(u, phi_x, phi_y, phi_x, phi_y, "Vertex vector fields scaled by lattice spacing")

if __name__ == "__main__":
    main()

