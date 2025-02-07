import numpy as np
import matplotlib.pyplot as plt
import scipy
from shapely.geometry import Point, Polygon


def bdy_fn_left(x, N, alpha, R):
    if -N - alpha*1/2 < x < R -N - alpha:
        return np.sqrt(R**2 - (x + alpha + N)**2)
    else:
        exit
    
def bdy_fn_right(x, N, alpha, R):
    if -N - alpha*1/2 < x < R - N:
        return np.sqrt(R**2 - (x + N)**2)
    else:
        exit
    

def make_grid(x_min, x_max, y_min, y_max, n_pts):
    x_pts = np.linspace(x_min, x_max, n_pts)
    y_pts = np.linspace(y_min,y_max,n_pts)
    grid_pts = np.zeros([n_pts**2, 3])
    for i in range(len(x_pts)):
        for j in range(len(x_pts)):
            grid_pts[i + j*n_pts, 0] =  x_pts[i]
            grid_pts[i + j*n_pts, 1] = y_pts[j]
    
    return grid_pts   
    
# For now, suppose the curve is defined by alpha and R.
def define_bdy_curve(alpha, R, N, x_min, x_max):
    """
    Creates points forming the outline of the 
    given curve. In this case, it's the intersection of two
    circular arcs, parametrized by alpha and R, with resolution
    defined by N.
    """
    x_vals = np.linspace(x_min, x_max, N)
    
    left_bdy = bdy_fn_left(x_vals, N, alpha, R)
    right_bdy = bdy_fn_right(x_vals, N, alpha, R)
    
    return left_bdy, right_bdy


def label_solid_pts(*bdy, grid_pts):
    x_pts, y_pts = grid_pts[:,0], grid_pts[:,1]
    return 



def main():
    
    alpha = 0.2
    R = 1
    N = 100
    n_pts = 200
    x_min, x_max = -2, 2
    y_min, y_max = 0, 4
    
    grid_pts = make_grid(x_min, x_max, y_min, y_max, n_pts)
    left_bdy, right_bdy = define_bdy_curve(alpha, R, N, x_min, x_max)
    
    return 
    
main()