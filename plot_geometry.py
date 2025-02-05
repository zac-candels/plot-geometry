import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import fsolve

plt.close('all')

class BdyPts:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    distances = []
    normals = []


alpha = 0.2
R = 1


x_vals = np.linspace(-1,1.1,100000)
n_pts = 10
x_min, x_max = -1, 1
y_min, y_max = 0, 2
x_pts = np.linspace(-1, 1.1, n_pts)
dx = x_pts[1] - x_pts[0]
y_pts = np.linspace(0,2,n_pts)
dy = y_pts[1] - y_pts[0]
grid_pts = np.zeros([n_pts**2, 3])

vels = np.array([
    [1, 0],
    [0, 1],
    [-1, 0],
    [0,-1],
    [np.sqrt(2)/2, np.sqrt(2)/2],
    [-np.sqrt(2)/2, np.sqrt(2)/2],
    [np.sqrt(2)/2, -np.sqrt(2)/2],
    [-np.sqrt(2)/2, -np.sqrt(2)/2]   ])


for i in range(len(x_pts)):
    for j in range(len(x_pts)):
        grid_pts[i + j*n_pts, 0] =  x_pts[i]
        grid_pts[i + j*n_pts, 1] = y_pts[j]

X, Y = np.meshgrid(x_pts, y_pts)

def f(t):
    return 0. + t*0 - np.sqrt(R**2 - (0.643 + t*1 + alpha)**2)
    

def bdy_curve_left(x, y, N):
    if -N - alpha*1/2 < x < R -N - alpha:
        return y - np.sqrt(R**2 - (x + alpha + N)**2)
    else:
        return

def bdy_curve_left_directional(t, x0, y0, v, N):
    v_x, v_y = v
    return bdy_curve_left(x0 + t*v_x, y0 + t*v_y, N)

def bdy_curve_right(x, y, N):
    if -N - alpha*1/2 < x < R - N:
        return y - np.sqrt(R**2 - (x + N)**2)
    else:
        return
    
def bdy_curve_right_directional(t, x0, y0, v, N):
    v_x, v_y = v
    return bdy_curve_right(x0 + t*v_x, y0 + t*v_y, N)

# Compute y-values
def fn_left(x, N):
    if -N - alpha*1/2 < x < R -N - alpha:
        return np.sqrt(R**2 - (x + alpha + N)**2)
    else:
        return 

def fn_right(x, N):
    if -N - alpha*1/2 < x < R - N:
        return np.sqrt(R**2 - (x + N)**2)
    else:
        return


N = 1

# Create boundary curves
y_l = np.zeros([len(x_vals), 5])
y_r = np.zeros([len(x_vals), 5])
for i in range(N):
    for j in range(len(x_vals)):
        y_l[j, i] = fn_left(x_vals[j], i)
        y_r[j, i] = fn_right(x_vals[j], i)
    
    
    
    
# Label solid nodes
fluid_pts = []
solid_pts = []
bdy_pts = []
for k in range(n_pts**2):
    if grid_pts[k,0] < np.sqrt(R**2 - grid_pts[k,1]**2):
        if grid_pts[k,0] > np.sqrt(R**2\
                                     - grid_pts[k,1]**2\
                                         - 2*grid_pts[k,0]*alpha - alpha**2):
            grid_pts[k,2] = 2
            solid_pts.append( [ grid_pts[k,0], grid_pts[k,1] ] )

solid_pts = np.asarray(solid_pts)
            
                
# Label boundary nodes
for i in range(len(grid_pts)):
    if grid_pts[i,2] == 2:
        continue
    else:
        p = np.array([ grid_pts[i,0], grid_pts[i,1] ])
        for j in range(len(solid_pts)):
            q = np.array([ solid_pts[j,0], solid_pts[j,1] ])
            if np.linalg.norm(p - q) < np.sqrt(2)*np.max([dx,dy]):
                grid_pts[i,2] = 1
                bdy_pts.append( [ grid_pts[i,0], grid_pts[i,1] ] )
                
bdy_pts = np.asarray(bdy_pts)


distances = []
t_vals = np.linspace(0,2,100) 
misc_vals = np.zeros(100)
# Loop through all boundary points, calculate
# distances to the wall and normal vectors to the bdy curves.
for i in range(len(bdy_pts)):
    for j in range(len(vels)):
        x0, y0 = bdy_pts[i, 0], bdy_pts[i, 1]
        if bdy_curve_left(x0, y0, N) == np.nan or bdy_curve_right(x0, y0, N) == np.nan:
            continue
        vel_vec = vels[j, :]
        #f = lambda t: bdy_curve_right_directional(t, x0, y0, vel_vec, N)
        

        sol = scipy.optimize.root_scalar(f, bracket=[0, R-alpha], method='brentq')
        




plt.figure()
for i in range(N):
    plt.plot(x_vals, y_l[:,i])
    plt.plot(x_vals, y_r[:,i] )
    plt.xlim([-(N-1)-alpha*1/2-0.1, R - (N-1) + 0.1])
    plt.axhline(0, color='black')
    #plt.scatter(X, Y, color="black", s=1)
    #plt.scatter(solid_pts[:,0],solid_pts[:,1], color="red", marker="o", s=1)
    #plt.scatter(bdy_pts[:,0],bdy_pts[:,1], color="green", marker="o", s=5)
    
    
    
    