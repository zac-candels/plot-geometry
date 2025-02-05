import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import fsolve

plt.close('all')

"""
class GridNode:
    def __init__(self, x, y, node_type):
        self.x = x
        self.y = y
        self.type = node_type

        
class BoundaryNode(GridNode):
    def __init__(self, x, y):
        super().__init__(x, y, "boundary")
    distances = [] #np.array(distances)
    normals = [] #np.array(normals)
"""

class BoundaryNode:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    distances = []
    normals = []


def make_grid(x_min, x_max, y_min, y_max, n_pts):
    
    x_pts = np.linspace(x_min, x_max, n_pts)
    y_pts = np.linspace(y_min,y_max,n_pts)
    grid_pts = np.zeros([n_pts**2, 3])
    for i in range(len(x_pts)):
        for j in range(len(x_pts)):
            grid_pts[i + j*n_pts, 0] =  x_pts[i]
            grid_pts[i + j*n_pts, 1] = y_pts[j]
    
    return grid_pts

    
# Compute y-values
def fn_left(x, N, alpha, R):
    if -N - alpha*1/2 < x < R -N - alpha:
        return np.sqrt(R**2 - (x + alpha + N)**2)
    else:
        return 

def fn_right(x, N, alpha, R):
    if -N - alpha*1/2 < x < R - N:
        return np.sqrt(R**2 - (x + N)**2)
    else:
        return



def label_solid_nodes(grid_pts, R, alpha):
    n_pts_tot = len(grid_pts[:,0])
    solid_pts = []
    for k in range(n_pts_tot):
        x, y = grid_pts[k,0], grid_pts[k,1]
        if R**2 - y**2 < 0:
            continue
        if R**2 - y**2 - 2*x*alpha - alpha**2 < 0:
            continue
        
        if x < np.sqrt(R**2 - y**2):
            if x > np.sqrt(R**2 - y**2 - 2*x*alpha - alpha**2):
                grid_pts[k,2] = 2
                print("solid pt found")
                solid_pts.append( [ x, y ] )
    
    solid_pts = np.asarray(solid_pts)
    return solid_pts
            
                
# Label boundary nodes
def label_bdy_nodes(grid_pts, solid_pts, R, alpha):
    x_pts = grid_pts[:,0]
    dx = x_pts[1] - x_pts[0]
    
    y_pts = np.unique(grid_pts[:,1])
    dy = y_pts[1] - y_pts[0]
    bdy_pts = []
    boundary_nodes = {}
    eps = 1e-5
    for i in range(len(grid_pts)):
        if grid_pts[i,2] == 2:
            continue
        else:
            p = np.array([ grid_pts[i,0], grid_pts[i,1] ])
            for j in range(len(solid_pts)):
                q = np.array([ solid_pts[j,0], solid_pts[j,1] ])
                if np.linalg.norm(p - q) <= np.sqrt(2)*np.max([dx,dy]) + eps:
                    grid_pts[i,2] = 1
                    x0, y0 = grid_pts[i,0], grid_pts[i,1]
                    bdy_pts.append( [ x0, y0 ] )
                    bdy_pt = BoundaryNode(x0, y0)
                    boundary_nodes[(x0, y0)] = bdy_pt
                    
    bdy_pts = np.asarray(bdy_pts)
    return boundary_nodes


# Loop through all boundary points, calculate
# distances to the wall and normal vectors to the bdy curves.
def distances_and_normals(boundary_nodes, R, alpha, vels, N):
    
    for (x0,y0), _ in boundary_nodes.items():
        distance_list = []
        for j in range(len(vels)):
            v_x, v_y = vels[j, :]
            coeffs = [v_x**2 + v_y**2, 2*(y0*v_x + x0*v_y + v_x*alpha),\
                      y0**2 - R**2 + x0**2 + alpha**2 + 2*x0*alpha]
            roots = np.roots(coeffs)
            if np.iscomplexobj(roots) == True:
                continue
            
            if len(roots) == 1 and roots > 0:
                distance_list.append(roots)
            elif len(roots) == 1 and roots < 0:
                continue
            elif len(roots) == 2 and np.all(roots > 0):
                distance_list.append(np.min(roots))
            elif len(roots) == 2 and np.all(roots < 0):
                continue
            elif len(roots) == 2 and np.any(roots < 0):
                distance_list.append( roots[roots> 0])
        
        #bdy_data.distances = distance_list
    
    
            
            
            #f = lambda t: bdy_curve_right_directional(t, x0, y0, vel_vec, N)
            
            


def visualize(x_vals, grid_pts, solid_pts, boundary_nodes, N, alpha, R):
    
    x, y = grid_pts[:,0], grid_pts[:,1]
    bdy_pts = np.asarray(list(boundary_nodes.keys()))
    y_l = np.zeros([len(x_vals), 5])
    y_r = np.zeros([len(x_vals), 5])
    for i in range(N):
        for j in range(len(x_vals)):
            y_l[j, i] = fn_left(x_vals[j], i, alpha, R)
            y_r[j, i] = fn_right(x_vals[j], i, alpha, R)
            
    plt.figure()
    for i in range(N):
        plt.plot(x_vals, y_l[:,i])
        plt.plot(x_vals, y_r[:,i] )
        plt.xlim([-(N-1)-alpha*1/2-0.1, R - (N-1) + 0.1])
        plt.axhline(0, color='black')
        plt.scatter(x, y, color="black", s=1)
        plt.scatter(solid_pts[:,0],solid_pts[:,1], color="red", marker="o", s=1)
        plt.scatter(bdy_pts[:,0],bdy_pts[:,1], color="green", marker="o", s=5)



def main():
    vels = np.array([
        [1, 0],
        [0, 1],
        [-1, 0],
        [0,-1],
        [np.sqrt(2)/2, np.sqrt(2)/2],
        [-np.sqrt(2)/2, np.sqrt(2)/2],
        [np.sqrt(2)/2, -np.sqrt(2)/2],
        [-np.sqrt(2)/2, -np.sqrt(2)/2]   ])
    
    N = 1 # Number of times boundary shape is repeated
    n_pts_per_direction = 10
    x_min, x_max = -1, 1
    y_min, y_max = 0, 2
    R = 1
    alpha = 0.2
    x_vals = np.linspace(x_min,x_max + 0.1,100000)
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_pts_per_direction)
    solid_pts = label_solid_nodes(grid_pts, R, alpha)
    boundary_nodes = label_bdy_nodes(grid_pts, solid_pts, R, alpha)
    
    distances_and_normals(boundary_nodes, R, alpha, vels, N)

    
    visualize(x_vals, grid_pts, solid_pts, boundary_nodes, N, alpha, R)
    
    
main()
    
    
    
    