# ----------------- Header -----------------------------------
# Author: Brigham Ostergaard
# Date: 8-15-2025
# Title: 1D_steady_conduction.py
# Description: Getting back to understanding how to do the CFD stuff. This is an 1D steady state conduction problem
#   Boundary conditions are set at each end of a rod. A source term is added as a constant (for now). We solve the temperature
#   profile of the rod using a*T = b. We solve for the a and b matricies, and use them to solve for T.
#-------------------------------------------------------------
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

# ------------------ Geometry/ Mesh setup ---------------------------------------
L = 5   # length of rod (m)
w = 0.25   # width of rod (m)
nx = 5   # nodes, element discretization
ny = 1

method = "a"            # This changes the node placements from type "a" to type "b" sets. I don't know if "b" works right, but it may be close
show_plots = True

dL = L/(nx-1)

# -------------- Physics Conditions/ Boundary Conditions ------------------------

# BC's
b1_method = "temp"      # define what type of condition
b2_method = "temp"

T0 = 80                # give each BC a value, Fixed temps in degrees C
TL = 20

# Physics
k1 = 205                   # conductive heat transfer coefficient, W/m (aluminum)
k2 = 109                   # brass
k2 = 0.1                   # plastic
transition = 0.5           # percentage of rod which is material 1

# ---------------------- Mesh Generation ----------------------------------------
if(method == "a"):
    points = np.linspace(dL/2,L-dL/2,nx-2)
    points = np.hstack([0,points])
    points = np.hstack([points,L])
elif(method == "b"):
    points = np.linspace(-dL/2,L+dL/2,nx)

# new_points = np.array([])
# for i in range(ny-1):
#     points = np.hstack([points,points])

# points = np.tile(points,ny)
# points = np.reshape(points,(ny,nx))

# ------------------------ Show mesh -------------------------------------------
lim = max([w,L])
if(show_plots):
    plt.figure()
    plt.plot([0,L,L,0,0],[w/2,w/2,-w/2,-w/2,w/2])       #gives outline of the bar (what we are simulating our 1D conduction)
    plt.scatter(points,np.zeros(np.size(points)),1,color='red')
    plt.ylim([-lim/10-lim,lim+lim/10])
    plt.xlim([-lim/10, lim+lim/10])



# define our "a" and "b" array/matrix
ap = np.zeros_like(points)
# np.fill_diagonal(ap,1)
aw = np.zeros_like(points)
ae = np.zeros_like(points)
b = np.transpose(np.zeros_like(points))

S = 0      # heat generation per m
nk1 = 0
nk2 = 0
k = k1

for i in range(nx):
    if(i > 0):
        if(method == "a"):
            if(i < nx-1):
                aw[i] = k*dL/w
            else:
                aw[i] = k*dL/(2*w)
        else:
            aw[i] = k*dL/w
    else:
        aw[i] = 0
    if(i < nx-1):
        if(method == "a"):
            if(i > 0):
                ae[i] = k*dL/w
            else:
                ae[i] = k*dL/(2*w)
        else:
            ae[i] = k*dL/w
    else:
        ae[i] = 0

    ap[i] = S/dL + (aw[i] + ae[i])      # ap - aw - ae - S = 0

# make "A" matrix  (A*T = b)
A = np.zeros([nx,nx])
for i in range(nx):
    if(i > 0):
        A[i,i-1] = aw[i]
    if(i < nx-1):
        A[i,i+1] = ae[i]
    A[i,i] = ap[i]

# make "b" vector, set with boundary conditions. (A*T = b)
b = np.zeros([nx,1])
b[0] = T0
b[nx-1] = TL

T = la.solve(A,b)

# print(T)

if(show_plots == True):
    plt.figure()
    plt.plot(points[0,:],T)
    plt.xlabel('Location (m)')
    plt.ylabel('Temp ($^\circ$C)')
    # plt.show()

    contour_map_y = np.linspace(w/2,-w/2,ny)

    # define x,y coordinates
    plot_plane_x = points[0,:]
    plot_plane_y = contour_map_y

    X,Y = np.meshgrid(plot_plane_x,plot_plane_y)
    reps = np.shape(plot_plane_y)
    T_plot = np.tile(np.array([T]),reps)
    # ny = np.size(contour_map_y)
    # nx = n
    T_plot = np.reshape(T_plot, (ny,nx))

    plt.figure()
    plt.contourf(X,Y,T_plot)
    plt.ylim([-lim/10-lim,lim+lim/10])
    plt.xlim([-lim/10, lim+lim/10])
    plt.show()