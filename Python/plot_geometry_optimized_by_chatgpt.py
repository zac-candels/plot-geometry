import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.interpolate import interp1d
import time

class BoundaryNode:
    def __init__(self, x, y, velocity_vecs):
        self.x = x
        self.y = y
        self.velocity_vecs = velocity_vecs
        self.distances = []
        self.normals = []

def make_grid(x_min, x_max, y_min, y_max, n_grid_pts):
    """Creates a uniform grid of points."""
    x_pts = np.linspace(x_min, x_max, n_grid_pts)
    y_pts = np.linspace(y_min, y_max, n_grid_pts)
    X, Y = np.meshgrid(x_pts, y_pts)
    return np.column_stack([X.ravel(), Y.ravel(), np.zeros(n_grid_pts**2)])

def define_bdy_curve(alpha, R, N_bdy_pts, x_min, x_max):
    """Defines the boundary curve as a sinusoidal function."""
    eps = 0.2
    x_vals = np.linspace(x_min - eps, x_max + eps, N_bdy_pts)
    y_vals = np.sin(5 * np.pi * x_vals) ** 2
    return np.column_stack([x_vals, y_vals])

def label_solid_pts(grid_pts, bdy_curve_pts):
    """Labels grid points that are below the boundary curve."""
    curve_interp = interp1d(bdy_curve_pts[:, 0], bdy_curve_pts[:, 1], kind='linear', fill_value="extrapolate")
    below_curve = grid_pts[:, 1] < curve_interp(grid_pts[:, 0])
    grid_pts[below_curve, 2] = 2  # Label solid points as '2'
    return grid_pts[below_curve, :2]  # Return only (x, y) of solid points

def label_bdy_points(grid_pts, solid_pts):
    """Finds boundary points by checking if a grid point is adjacent to a solid point."""
    kdtree = KDTree(solid_pts)
    dists, _ = kdtree.query(grid_pts[:, :2], k=1)
    
    is_boundary = (dists <= np.sqrt(2) * (grid_pts[1, 0] - grid_pts[0, 0])) & (grid_pts[:, 2] != 2)
    grid_pts[is_boundary, 2] = 1  # Label boundary points as '1'

    boundary_nodes = []
    for x, y in grid_pts[is_boundary, :2]:
        neighbors = solid_pts[np.linalg.norm(solid_pts - [x, y], axis=1) <= np.sqrt(2) * (grid_pts[1, 0] - grid_pts[0, 0])]
        velocity_vecs = (neighbors - [x, y]) / np.linalg.norm(neighbors - [x, y], axis=1)[:, None]
        boundary_nodes.append(BoundaryNode(x, y, velocity_vecs.tolist()))
    
    return boundary_nodes

def distances_and_normals(boundary_nodes, bdy_curve_pts):
    """Computes distances and normals to the boundary."""
    kdtree = KDTree(bdy_curve_pts)
    
    for bdy_node in boundary_nodes:
        x_b = np.array([bdy_node.x, bdy_node.y])
        for unit_vel_vec in bdy_node.velocity_vecs:
            _, idx = kdtree.query(x_b + unit_vel_vec)  # Nearest boundary point in direction of velocity
            pt2 = bdy_curve_pts[idx]
            normal = np.array([-unit_vel_vec[1], unit_vel_vec[0]])  # Perpendicular to velocity
            
            delta = np.linalg.norm(pt2 - x_b)  # Approximate distance
            bdy_node.distances.append(delta)
            bdy_node.normals.append(normal / np.linalg.norm(normal))
    
    return boundary_nodes

def visualize(grid_pts, bdy_curve_pts, solid_pts, bdy_nodes):
    """Visualizes the grid, boundary, and classified points."""
    bdy_pts = np.array([[node.x, node.y] for node in bdy_nodes])

    plt.figure(figsize=(8, 6))
    plt.scatter(grid_pts[:, 0], grid_pts[:, 1], s=0.1, color='gray', label="Grid Points")
    plt.plot(bdy_curve_pts[:, 0], bdy_curve_pts[:, 1], color='black', label="Boundary Curve")
    plt.scatter(solid_pts[:, 0], solid_pts[:, 1], marker="x", color="red", label="Solid Points")
    plt.scatter(bdy_pts[:, 0], bdy_pts[:, 1], color="green", label="Boundary Points")
    plt.legend()
    plt.show()

def main():
    """Main function to run the boundary classification process."""
    start_time = time.time()
    
    # Parameters
    alpha, R = 0.2, 1
    N_bdy_pts = 3200
    n_grid_pts = 100
    x_min, x_max, y_min, y_max = -2, 2, 0, 4
    
    # Step 1: Create grid
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts)
    
    # Step 2: Define boundary curve
    bdy_curve_pts = define_bdy_curve(alpha, R, N_bdy_pts, x_min, x_max)
    
    # Step 3: Label solid points
    solid_pts = label_solid_pts(grid_pts, bdy_curve_pts)
    
    # Step 4: Label boundary points
    bdy_nodes = label_bdy_points(grid_pts, solid_pts)
    
    # Step 5: Compute distances and normals
    bdy_nodes = distances_and_normals(bdy_nodes, bdy_curve_pts)
    
    # Step 6: Visualize results
    visualize(grid_pts, bdy_curve_pts, solid_pts, bdy_nodes)
    
    print("Run time:", time.time() - start_time)

main()
