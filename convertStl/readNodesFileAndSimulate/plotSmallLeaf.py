import numpy as np
import math
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

# --- Parameters ---
filename = '/home/zcandels/plot-geometry/convertStl/readNodesFileAndSimulate/smallLeaf.txt'

TX = 74; TY = 29; TZ = 43

PX = 10
PY = 1
PZ = 1

#smallLeafFile = open("smallLeaf.txt", "w")

# --- Load data ---
v = np.fromfile(filename, dtype=np.int32, sep=' ')
# Reshape: MATLAB reshapes column-major (Fortran order)
G = v.reshape((TX, TY, TZ), order='C')

plt.figure()
plt.imshow(G[30,:,:].T, origin='lower')
plt.title("x slice")
plt.show()

plt.figure()
plt.imshow(G[:,10,:].T, origin='lower')
plt.title("y slice")
plt.show()

plt.figure()
plt.imshow(G[:,:,20].T, origin='lower')
plt.title("z slice")
plt.show()

# --- Write partitioned data files ---
os.makedirs('data', exist_ok=True)

ranks = []
starts = []
lengths = []

for a in range(1, PX * PY * PZ + 1):
    rankz = math.ceil(a / (PX * PY))
    ranky = math.ceil((a - (rankz - 1) * PX * PY) / PX)
    rankx = a - (rankz - 1) * PX * PY - (ranky - 1) * PX

    ranks.append((rankx, ranky, rankz))

    # X partition
    if rankx <= TX % PX:
        lengthx = math.ceil(TX / PX)
        startx  = (rankx - 1) * lengthx
    else:
        lengthx = math.floor(TX / PX)
        startx  = (rankx - 1) * lengthx + TX % PX

    # Y partition
    if ranky <= TY % PY:
        lengthy = math.ceil(TY / PY)
        starty  = (ranky - 1) * lengthy
    else:
        lengthy = math.floor(TY / PY)
        starty  = (ranky - 1) * lengthy + TY % PY

    # Z partition
    if rankz <= TZ % PZ:
        lengthz = math.ceil(TZ / PZ)
        startz  = (rankz - 1) * lengthz
    else:
        lengthz = math.floor(TZ / PZ)
        startz  = (rankz - 1) * lengthz + TZ % PZ

    starts.append((startx, starty, startz))
    lengths.append((lengthx, lengthy, lengthz))

    NX = lengthx + 2
    NY = lengthy + 2
    NZ = lengthz + 2

    data = np.zeros((NX, NY, NZ), dtype=np.int32)

    for i in range(1, NX + 1):
        for j in range(1, NY + 1):
            for k in range(1, NZ + 1):
                # Periodic boundary in X
                id_ = startx + i - 1
                if id_ == 0:
                    id_ = TX
                if id_ == TX + 1:
                    id_ = 1

                # Periodic boundary in Y
                jd = starty + j - 1
                if jd == 0:
                    jd = TY
                if jd == TY + 1:
                    jd = 1

                # Periodic boundary in Z
                kd = startz + k - 1
                if kd == 0:
                    kd = TZ
                if kd == TZ + 1:
                    kd = 1

                # Convert to 0-based indexing for numpy
                data[i - 1, j - 1, k - 1] = G[id_ - 1, jd - 1, kd - 1]

    mpirank_str = f'{a - 1:04d}'
    out_filename = f'./data/data{mpirank_str}.dat'
    with open(out_filename, 'w') as fid:
        for i in range(NX):
            for j in range(NY):
                for k in range(NZ):
                    fid.write(f'{int(data[i, j, k])}\n')