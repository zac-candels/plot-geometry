import numpy as np
import math
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

# --- Parameters ---
filename = '/home/zcandels/plot-geometry/convertStl/readNodesFileAndSimulate/nodes.out'

TX = 764; TY = 29; TZ = 82

PX = 10
PY = 1
PZ = 1

smallLeafFile = open("smallLeaf.txt", "w")
largeLeafFile = open("largeLeaf.txt", "w")

# --- Load data ---
v = np.fromfile(filename, dtype=np.float32, sep=' ')

assert v.size == TX * TY * TZ, f'数据量不匹配：期望 {TX*TY*TZ}，实际 {v.size}'

# Reshape: MATLAB reshapes column-major (Fortran order)
tmp = v.reshape((TZ, TY, TX), order='F')   # tmp[k, j, i]

A = tmp.copy()  # A[k, j, i], shape: (TZ, TY, TX)

TZ = 43  # Truncate Z dimension for loop below

print(A.shape)
print(f'{A.min():.6f}  {A.max():.6f}')

# --- Build G: reverse X axis, swap k and i ---
# MATLAB: G(TX-i+1, j, k) = A(k, j, i)  for i=1..TX, j=1..TY, k=1..TZ
# In 0-based: G[TX-1-i, j, k] = A[k, j, i]
# This is equivalent to: take A[0:TZ, :, :] and transpose/flip axes

G = np.zeros((TX, TY, TZ), dtype=np.float32)
for i in range(TX):
    for j in range(TY):
        for k in range(TZ):
            G[TX - 1 - i, j, k] = A[k, j, i]

# --- Remap solid values ---
# Replace 1 -> 2 (solid with contact angle 30)
G[G == 1] = 2

# --- Liquid initialisation ---
# Replace -1 -> 1
G[G == -1] = 1

print("------------------")
print(G.shape)
print(np.unique(G))
for i in range(532, 606):
    for j in range(TY):
        for k in range(TZ):
            globalIndex = (i-1) + (j-1)*TX + (k-1)*TX*TY
            smallLeafFile.write(str(int(G[i,j,k])))
            smallLeafFile.write("\n")
            
ctr = 0
for i in range(600, 660):
    for j in range(TY):
        for k in range(TZ):

            if ( i  > 630 ) and (i <= 660):
                if ( j >= 10 ) and (j < 21):
                    if (k >=0) and (k <= 22):
                        if G[i,j,k] == 1:
                            largeLeafFile.write(str(2))
                            largeLeafFile.write("\n")
                            ctr +=1 
                            continue
            if (i > 652) and (i <= 660):
                if (j>=15) and (j <= TY):
                    if (k >= 0) and (k <= 15):
                        if G[i,j,k] == 1:
                            largeLeafFile.write(str(2))
                            largeLeafFile.write("\n")
                            ctr +=1 
                            continue
            
            largeLeafFile.write(str(int(G[i,j,k])))
            largeLeafFile.write("\n")
            ctr+=1
print("ctr = ", ctr)

smallLeafFile.close()
# --- Extend domain in X by Nx_extra ---
Nx_extra = 100
G_ext = np.ones((TX + Nx_extra, TY, TZ), dtype=np.float32)

# MATLAB: G_ext(Nx_extra : Nx_extra-1+TX, :, :) = G
# MATLAB 1-based inclusive [Nx_extra, Nx_extra-1+TX] -> 0-based [Nx_extra-1, Nx_extra-1+TX)
G_ext[Nx_extra - 1 : Nx_extra - 1 + TX, :, :] = G

TX_old = TX
TX = TX + Nx_extra

G = G_ext

# --- Set a block of cells to solid (2) ---
# MATLAB: i=1:100, j=1:10  -> 0-based: i=0:100, j=0:10
G[0:100, 0:10, 0:TZ] = 2

# --- Visualisation (isosurface equivalent) ---
try:
    verts, faces, _, _ = measure.marching_cubes(G, level=0.5)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(verts[faces], alpha=0.3)
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    ax.set_xlim(0, G.shape[0])
    ax.set_ylim(0, G.shape[1])
    ax.set_zlim(0, G.shape[2])
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()
except Exception as e:
    print(f'Visualisation skipped: {e}')

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