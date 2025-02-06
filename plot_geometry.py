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

    
# Needed for visualization
def bdy_fn_left(x, N, alpha, R):
    if -N - alpha*1/2 < x < R -N - alpha:
        return np.sqrt(R**2 - (x + alpha + N)**2)
    else:
        return 

def deriv_bdy_fn_left(x, y, N, alpha, R):
    return -(x + alpha)/y

def bdy_fn_right(x, N, alpha, R):
    if -N - alpha*1/2 < x < R - N:
        return np.sqrt(R**2 - (x + N)**2)
    else:
        return

def deriv_bdy_fn_right(x, y, N):
    return -x/y



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
                #print("solid pt found")
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
    
    for (x0,y0), bdy_node in boundary_nodes.items():
        if np.abs(y0 - 0) < 1e-4:
            continue
        distance_list = []
        normals_list = []
        for j in range(len(vels)):
            v_x, v_y = vels[j, :]
            # roots1 corresponds to left fn
            coeffs1 = [v_x**2 + v_y**2, 2*(y0*v_x + x0*v_y + v_x*alpha),\
                      y0**2 - R**2 + x0**2 + alpha**2 + 2*x0*alpha]
            roots1 = np.roots(coeffs1)
            
            # roots2 corresponds to right fn
            coeffs2 = [v_x**2 + v_y**2, 2*(y0*v_x + x0*v_y),\
                      y0**2 - R**2 + x0**2]
            roots2 = np.roots(coeffs2)
            
            if np.any(roots1 < 0) == True and np.any(roots1 > 0) == True:
                roots1 = roots1[roots1 > 0]
                
            if np.any(roots2 < 0) == True and np.any(roots2 > 0) == True:
                roots2 = roots2[roots2 > 0]
                
            if np.all(roots1 > 0):
                roots1 = np.min(roots1)
                
            if np.all(roots2 > 0):
                roots2 = np.min(roots2)
                
            if np.all(roots1 < 0) and np.all(roots2 < 0):
                continue

            if np.iscomplexobj(roots1) == True\
                and np.iscomplexobj(roots2) == True:
                continue

            if np.iscomplexobj(roots1) == True\
                and np.iscomplexobj(roots2) == False\
                    and np.all(roots2< 0):
                        continue

            if np.iscomplexobj(roots2) == True\
                and np.iscomplexobj(roots1) == False\
                    and np.all(roots1< 0):
                        continue

            if np.iscomplexobj(roots1) == True\
                and np.iscomplexobj(roots2) == False\
                    and np.all(roots2> 0):
                        print("dist found")
                        which_fn = "Right"
                        dist = roots2
                        distance_list.append(dist)
                        
            if np.iscomplexobj(roots2) == True\
                and np.iscomplexobj(roots1) == False\
                    and np.all(roots1> 0):
                        print("dist found")
                        which_fn = "Left"
                        dist = roots1
                        distance_list.append(dist)
                        
            if np.all(roots1 > 0) == True and np.all(roots2 > 0) == True:
                print("dist found")
                dist = np.min( [roots1, roots2] )
                if roots1 < roots2:
                    which_fn = "Left"
                else:
                    which_fn = "Right"
                distance_list.append( dist )
            
            if which_fn == "Right": # ie right bdy function
                slope_at_bdy_curve = deriv_bdy_fn_right(x0, y0, N)
                slope_of_normal = -1/slope_at_bdy_curve
                normal_vec = np.array( [1, slope_of_normal] )
                normal_vec = normal_vec/np.linalg.norm(normal_vec)
                normals_list.append(normal_vec)
                
            if which_fn == "Left": # ie right bdy function
                slope_at_bdy_curve = deriv_bdy_fn_left(x0, y0, N, alpha, R)
                slope_of_normal = -1/slope_at_bdy_curve
                normal_vec = np.array( [1, slope_of_normal] )
                normal_vec = normal_vec/np.linalg.norm(normal_vec)
                normals_list.append(normal_vec)
        bdy_node.distances = distance_list
        bdy_node.normals = normals_list
            
            
            
        
        #bdy_data.distances = distance_list


def visualize(x_vals, grid_pts, solid_pts, boundary_nodes, N, alpha, R):
    
    x, y = grid_pts[:,0], grid_pts[:,1]
    bdy_pts = np.asarray(list(boundary_nodes.keys()))
    y_l = np.zeros([len(x_vals), 5])
    y_r = np.zeros([len(x_vals), 5])
    for i in range(N):
        for j in range(len(x_vals)):
            y_l[j, i] = bdy_fn_left(x_vals[j], i, alpha, R)
            y_r[j, i] = bdy_fn_right(x_vals[j], i, alpha, R)
            
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
    
    
    
    