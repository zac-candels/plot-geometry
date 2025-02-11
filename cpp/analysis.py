import numpy as np 
import matplotlib.pyplot as plt 

path = "/home/zcandels/geom/cpp/"
grid_pts = np.loadtxt(path+"grid_pts.dat", delimiter=",")
bdy_curve_pts = np.loadtxt(path+"bdy_curve_pts.dat", delimiter=",")
solid_pts = np.loadtxt(path+"solid_pts.dat", delimiter=",")
bdy_pts = np.loadtxt(path+"boundary_pts.dat", delimiter=",")

plt.figure()
plt.scatter(grid_pts[:,0], grid_pts[:,1],s=0.1)
plt.plot(bdy_curve_pts[:,0], bdy_curve_pts[:,1])
plt.scatter(solid_pts[:,0], solid_pts[:,1], marker="x", color="red")
plt.scatter(bdy_pts[:,0], bdy_pts[:,1], color="green")
plt.show()