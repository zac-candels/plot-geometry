import numpy as np
import matplotlib.pyplot as plt
import scipy
from shapely.geometry import Point, Polygon


def bdy_fn_left(x, N_repeats, alpha, R):
    x = np.asarray(x)
    x = x[ (x < R - N_repeats - alpha) & (x > -N_repeats - alpha*1/2)]
    return np.sqrt(R**2 - (x + alpha + N_repeats)**2)
    
        
    
def bdy_fn_right(x, N_repeats, alpha, R):
    x = np.asarray(x)
    
    x = x[(x > - N_repeats - alpha*1/2) & (x < R - N_repeats)]
    return np.sqrt(R**2 - (x + N_repeats)**2)
    

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


def label_solid_pts(*bdy, grid_pts, left_bdy_pts, right_bdy_pts):
    x_pts, y_pts = grid_pts[:,0], grid_pts[:,1]
    x_min, x_max = x_pts[0], x_pts[-1]
    y_min, y_max = y_pts[0], y_pts[-1]
    
    
    
    return 



def main():
    
    alpha = 0.2
    R = 1
    N_bdy_pts = 100
    N_repeats = 1
    n_grid_pts = 200
    x_min, x_max = -2, 2
    y_min, y_max = 0, 4
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts)
    
    left_bdy_pts, right_bdy_pts = define_bdy_curve(alpha, R,
                                           N_bdy_pts, N_repeats, x_min, x_max)
    
    return 
    
main()