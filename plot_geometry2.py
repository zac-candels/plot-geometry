import numpy as np
import matplotlib.pyplot as plt
import scipy
from shapely.geometry import Point, Polygon

plt.close('all')


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
    for point in mesh_pts:
        if bdy_shape_polygon.contains(Point(point)):
            pts_x.append(point[0])
            pts_y.append(point[1])
            
            solid_points.append(point)
    solid_points = np.transpose( np.asarray( [ pts_x, pts_y ] ) )
    return solid_points


def visualize(grid_pts, left_bdy_pts, right_bdy_pts, solid_pts):
    plt.figure()
    plt.scatter(grid_pts[:,0], grid_pts[:,1])
    plt.plot(left_bdy_pts[:,0], left_bdy_pts[:,1])
    plt.plot(right_bdy_pts[:,0], right_bdy_pts[:,1])
    plt.scatter(solid_pts[:,0], solid_pts[:,1], color="red")
    

def main():
    
    alpha = 0.2
    R = 1
    N_bdy_pts = 100
    N_repeats = 1
    n_grid_pts = 10
    x_min, x_max = -2, 2
    y_min, y_max = 0, 4
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts)
    
    left_bdy_pts, right_bdy_pts = define_bdy_curve(alpha, R,
                                           N_bdy_pts, N_repeats, x_min, x_max)
    
    solid_points = label_solid_pts(grid_pts, left_bdy_pts, right_bdy_pts)
    
    visualize(grid_pts, left_bdy_pts, right_bdy_pts, solid_points)
    
    return 
    
main()