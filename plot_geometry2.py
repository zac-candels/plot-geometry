import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
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
    distances = []
    normals = []


def bdy_fn_left(x, N_repeats, alpha, R):
    x = np.asarray(x)
    x = x[ (x < R - N_repeats - alpha) & (x > -N_repeats - alpha*1/2)]
    y = np.sqrt(R**2 - (x + alpha + N_repeats)**2)
    left_bdy_pts = np.transpose(np.array([x, y]))
    return left_bdy_pts
    
        
    
def bdy_fn_right(x, N_repeats, alpha, R):
    x = np.asarray(x)
    
    x = x[(x > - N_repeats - alpha*1/2) & (x < R - N_repeats)]
    y = np.sqrt(R**2 - (x + N_repeats)**2)
    right_bdy_pts = np.transpose(np.array([x, y]))
    return right_bdy_pts
    

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
    x_vals = np.linspace(x_min, x_max, N_bdy_pts)
    
    left_bdy = bdy_fn_left(x_vals, N_repeats, alpha, R)
    right_bdy = bdy_fn_right(x_vals, N_repeats, alpha, R)
    
    return left_bdy, right_bdy


def label_solid_pts(grid_pts, left_bdy_pts, right_bdy_pts):
    
    mesh_pts = [tuple(row) for row in grid_pts[:,0:2]]
    right_bdy_pts = [tuple(row) for row in right_bdy_pts]
    left_bdy_pts = [tuple(row) for row in left_bdy_pts]
    
    bdy_shape_pts = left_bdy_pts + right_bdy_pts
    
    bdy_shape_polygon = Polygon(bdy_shape_pts)

    # Label points inside the crescent moon
    solid_points = []
    pts_x, pts_y = [], []
    for idx, point in enumerate(mesh_pts):
        if bdy_shape_polygon.contains(Point(point)):
            pts_x.append(point[0])
            pts_y.append(point[1])
            
            solid_points.append(point)
            grid_pts[idx,2] = 2
    solid_points = np.transpose( np.asarray( [ pts_x, pts_y ] ) )
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


def distances_and_normals(bdy_nodes, left_bdy_pts, right_bdy_pts):
    
    
    
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
        
        for idx in idxs:
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
        
    pt_set = np.concatenate( [left_bdy_pts, right_bdy_pts] )
    
    
    for location, node_info in bdy_nodes.items():
        x_b = np.asarray( location )
        for i in range(len(node_info.velocity_vecs)):
            unit_vel_vec = node_info.velocity_vecs[i]
            pt1, pt2, d = find_closest_points(pt_set, x_b, unit_vel_vec)  
            v = pt2 - pt1
            
            intersection_mat = np.array([ [v[0], -unit_vel_vec[0]],
                                         [v[1], -unit_vel_vec[1]] ])
            intersection_rhs_vec = np.array([ x_b[0] - pt1[0],
                                              x_b[1] - pt1[1] ])
            
            intersection_distances = np.linalg.solve(intersection_mat,
                                                     intersection_rhs_vec)
            
            delta = intersection_distances[1]
            if delta > 1:
                continue
            else:
                node_info.distances.append(delta)
            
            
            
        





def visualize(grid_pts, left_bdy_pts, right_bdy_pts, solid_pts, bdy_nodes):
    bdy_pts = np.asarray(list(bdy_nodes.keys()))
    
    plt.figure()
    plt.scatter(grid_pts[:,0], grid_pts[:,1])
    plt.plot(left_bdy_pts[:,0], left_bdy_pts[:,1])
    plt.plot(right_bdy_pts[:,0], right_bdy_pts[:,1])
    plt.scatter(solid_pts[:,0], solid_pts[:,1], color="red")
    plt.scatter(bdy_pts[:,0], bdy_pts[:,1], color="green")
    

def main():
    
    alpha = 0.2
    R = 1
    N_bdy_pts = 400
    N_repeats = 1
    n_grid_pts = 10
    x_min, x_max = -2, 2
    y_min, y_max = 0, 4
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts)
    
    left_bdy_curve_pts, right_bdy_curve_pts = define_bdy_curve(alpha, R,
                                           N_bdy_pts, N_repeats, x_min, x_max)
    
    solid_pts = label_solid_pts(grid_pts, left_bdy_curve_pts,
                                right_bdy_curve_pts)
    
    bdy_nodes = label_bdy_points(grid_pts, solid_pts)
    
    distances_and_normals(bdy_nodes, left_bdy_curve_pts, right_bdy_curve_pts)
    
    visualize(grid_pts, left_bdy_curve_pts,
              right_bdy_curve_pts, solid_pts, bdy_nodes)
    
    return 
    
main()