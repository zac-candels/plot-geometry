import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.interpolate import interp1d
from shapely.geometry import Point, Polygon

# This program does several things.
# 1. Contains a function to define the geometry of the system 
# and create points for the boundary curve(s)
# 2. Function to label the solid points.
# 3. Function to label boundary points
# 4. Function to compute distances and normals from boundary points
# to the boundary curve. 


# Note, to determine the position of a grid point in the c++ code,
# use the function ComputeXYZ in the service.hh header file.

#plt.close('all')

class BoundaryNode:
    def __init__(self, x, y, velocity_vecs):
        self.x = x
        self.y = y
        self.velocity_vecs = velocity_vecs
        self.distances = []
        self.normals = []
    

def make_grid(x_min, x_max, y_min, y_max, n_grid_pts):
    x_pts = np.linspace(x_min, x_max, n_grid_pts)
    y_pts = np.linspace(y_min,y_max,n_grid_pts)
    grid_pts = np.zeros([n_grid_pts**2, 3])
    for i in range(len(x_pts)):
        for j in range(len(x_pts)):
            grid_pts[i + j*n_grid_pts, 0] =  x_pts[i]
            grid_pts[i + j*n_grid_pts, 1] = y_pts[j]
    
    return grid_pts   
    
# For now, suppose the curve is defined by alpha and R.
def define_bdy_curve(alpha, R, N_bdy_pts, N_repeats, x_min, x_max):
    """
    Creates points forming the outline of the 
    given curve. In this case, it's the intersection of two
    circular arcs, parametrized by alpha and R, with resolution
    defined by N.
    """
    eps = 0.2
    x_vals = np.linspace(x_min - eps, x_max + eps, N_bdy_pts)
    
    y_vals = np.sin(5*np.pi*x_vals)*np.sin(5*np.pi*x_vals)
    
    bdy_curve_pts = np.transpose( [x_vals, y_vals] )
    
    return bdy_curve_pts, x_vals


def label_solid_pts(grid_pts, bdy_curve_pts):
    
    mesh_pts = [tuple(row) for row in grid_pts[:,0:2]]

    # Label points below bdy curve
    
    curve_interp = interp1d(bdy_curve_pts[:,0], bdy_curve_pts[:,1],
                            kind='linear', fill_value="extrapolate")
    
    solid_pts_x, solid_pts_y = [], []
    # x_mesh_pts, y_mesh_pts = grid_pts[:,0], grid_pts[:,1]
    # for i in range(len(x_mesh_pts)):
    #     for j in range(len(y_mesh_pts)):
    #         y_curve = curve_interp(x_mesh_pts[i])
    #         if y_mesh_pts[j] < y_curve:
    #             solid_pts_x.append(x_mesh_pts[i])
    #             solid_pts_y.append(y_mesh_pts[j])
                
    for idx, point in enumerate(mesh_pts):
        x_mesh_pt, y_mesh_pt = point[0], point[1]
        y_curve = curve_interp(point[0])
        if y_mesh_pt < y_curve:
            solid_pts_x.append(x_mesh_pt)
            solid_pts_y.append(y_mesh_pt)
            grid_pts[idx,2] = 2
            
    solid_points = np.transpose( np.asarray( [ solid_pts_x, solid_pts_y ] ) )
    return solid_points


def label_bdy_points(grid_pts, solid_pts):
    x_pts = grid_pts[:,0]
    dx = x_pts[1] - x_pts[0]
    
    y_pts = np.unique(grid_pts[:,1])
    dy = y_pts[1] - y_pts[0]
    bdy_pts = []
    boundary_nodes = {}
    eps = 1e-5
    for i in range(len(grid_pts)):
        velocity_vecs = []
        if grid_pts[i,2] == 2:
            continue
        else:
            P = np.array([ grid_pts[i,0], grid_pts[i,1] ])
            for j in range(len(solid_pts)):
                Q = np.array([ solid_pts[j,0], solid_pts[j,1] ])
                if np.linalg.norm(Q-P) <= np.sqrt(2)*np.max([dx,dy]) + eps:
                    velocity_vecs.append( (Q - P)/np.linalg.norm(Q-P) )
                    grid_pts[i,2] = 1
                    x0, y0 = grid_pts[i,0], grid_pts[i,1]
                    bdy_pts.append( [ x0, y0 ] )
                    for k in range(len(solid_pts)):
                        if k == j:
                            continue
                        Q = np.array([ solid_pts[k,0], solid_pts[k,1] ])
                        if np.linalg.norm(Q - P)\
                            <= np.sqrt(2)*np.max([dx,dy]) + eps:
                                
                           velocity_vecs.append( (Q - P)/np.linalg.norm(Q - P) ) 
                            
                    bdy_pt = BoundaryNode(x0, y0, velocity_vecs)
                    boundary_nodes[(x0, y0)] = bdy_pt
                    break
                    
    return boundary_nodes


def distances_and_normals(bdy_nodes, bdy_curve_pts):
    
    
    
    def signed_distance(p, x0, u):
        """Compute the signed distance of point p from the line x0 + t*u."""
        v_perp = np.array([-u[1], u[0]])  # Perpendicular unit vector
        return np.dot(p - x0, v_perp)
    
    
    
    def find_closest_points(pt_set, x0, unit_vec, k=10):
        debug_distances = []
        """Find the closest points in pt_set on either side of the line x0 + t*u."""
        # Build KD-tree for fast nearest neighbor search
        tree = KDTree(pt_set)
        
        # Find k nearest neighbors to x0
        _, idxs = tree.query(x0, k=k)
        
        closest_pos = None
        closest_neg = None
        min_pos_dist = float('inf')
        min_neg_dist = float('inf')
        
        for idx in range(len(pt_set)): #idxs:
            p = pt_set[idx]
            d = signed_distance(p, x0, unit_vec)
            
            debug_distances.append(d)
            
            if d > 0 and d < min_pos_dist:
                closest_pos = p
                min_pos_dist = d
            elif d < 0 and abs(d) < min_neg_dist:
                closest_neg = p
                min_neg_dist = abs(d)
        
        return closest_neg, closest_pos, debug_distances
        
    bdy_curve_pt_set = bdy_curve_pts
    
    
    for location, node_info in bdy_nodes.items():
        x_b = np.asarray( location )
        for i in range(len(node_info.velocity_vecs)):
            unit_vel_vec = node_info.velocity_vecs[i]
            pt1, pt2, d = find_closest_points(bdy_curve_pt_set,
                                              x_b, unit_vel_vec)  
            v = pt2 - pt1
            v = v/np.linalg.norm(v)
            
            intersection_mat = np.array([ [v[0], -unit_vel_vec[0]],
                                         [v[1], -unit_vel_vec[1]] ])
            intersection_rhs_vec = np.array([ x_b[0] - pt1[0],
                                              x_b[1] - pt1[1] ])
            
            intersection_distances = np.linalg.solve(intersection_mat,
                                                     intersection_rhs_vec)
            
            delta = intersection_distances[1]
            if delta > 1:
                print(location)
                continue
            else:
                node_info.distances.append(delta)
                normal = np.array( [ -1, -v[0]/v[1] ] )
                normal = normal/np.linalg.norm(normal)
                node_info.normals.append( np.array( [ 1, -v[0]/v[1] ] ) )
    
    return bdy_nodes
            



def visualize(grid_pts, bdy_curve_pts, solid_pts, bdy_nodes):
    bdy_pts = np.asarray(list(bdy_nodes.keys()))
    
    plt.figure()
    plt.scatter(grid_pts[:,0], grid_pts[:,1])
    plt.plot(bdy_curve_pts[:,0], bdy_curve_pts[:,1])
    plt.scatter(solid_pts[:,0], solid_pts[:,1], marker="x", color="red")
    plt.scatter(bdy_pts[:,0], bdy_pts[:,1], color="green")
    

def main():
    
    alpha = 0.2
    R = 1
    N_bdy_pts = 3200
    N_repeats = 1
    n_grid_pts = 10
    x_min, x_max = -2, 2
    y_min, y_max = 0, 4
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts)
    
    bdy_curve_pts, x_vals = define_bdy_curve(alpha, R,
                                           N_bdy_pts, N_repeats, x_min, x_max)
    
    solid_pts = label_solid_pts(grid_pts, bdy_curve_pts)
    
    bdy_nodes = label_bdy_points(grid_pts, solid_pts)
    
    bdy_nodes = distances_and_normals(bdy_nodes, bdy_curve_pts)
    
    visualize(grid_pts, bdy_curve_pts, solid_pts, bdy_nodes)
    
    return 
    
main()