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
L = 0.5   # length of rod (m)
w = 0.1   # diamter of rod (m) (used in CV sizing too, so we leave as w)
nx = 20   # nodes, element discretization
ny = 1

method = "a"            # This changes the node placements from type "a" to type "b" sets. I don't know if "b" works right, but it may be close
show_plots = True

if(nx > 3):
    dL = L/(nx-1)
# elif(nx == 2):
#     dL = L/3
else:
    dL = L/2
if(ny > 1):
    dw = w/(ny-1)
else:
    dw = w/2

# -------------- Physics Conditions/ Boundary Conditions ------------------------

# BC's
b1_method = "temp"      # define what type of condition
b2_method = "insulated"

T0 = 200                # give each BC a value, Fixed temps in degrees C
TL = 500

if(b1_method == "flux"):
    qL = 10                # prescribed flux at ends (only for "flux" bc)
elif(b1_method == "insulated"):
    qL = 0

if(b2_method == "flux"):
    qR = 10
elif(b2_method == "insulated"):
    qR = 0

# Physics
k1 = 205                   # conductive heat transfer coefficient, W/m (aluminum)
k1 = 1000                   # copper, or something. Example 4.1 from the book.
k2 = 109                   # brass
k2 = 0.1                   # plastic
transition = 0.5           # percentage of rod which is material 1

# Pipe Geometry
d = w
A_c = (np.pi/4)*(d**2)                  #A_c of the beam

# ---------------------- Mesh Generation ----------------------------------------
# create nodes and fill with x positions.
if(method == "a" and nx > 3):
    points_x = np.linspace(dL/2,L-dL/2,nx-2)
    points_x = np.hstack([0,points_x])
    points_x = np.hstack([points_x,L])
elif(method == "a" and nx == 3):
    points_x = np.array([0, L/2, L])
elif(method == "b"):
    points_x = np.linspace(-dL/2,L+dL/2,nx)

# if(method == "a" and ny > 3):
#     points_y = np.linspace(dw/2,w-dw/2,ny-2)
#     points_y = np.hstack([0,points_y])
#     points_y = np.hstack([points_y,w])
# elif(method == "a" and ny == 1):
#     points_y = np.array([dw])
# elif(method == "a" and ny == 2):
#     points_y = np.array([np.array([dw/2, w - dw/2])])
# elif(method == "a" and ny == 3):
#     points_y = np.array([np.array([0, w/2, w])])
# elif(method == "b" or ny <= 2):
#     points_y = np.linspace(-dw/2,w+dw/2,ny)

points_y = np.ones(ny)*w/2

points_x,points_y = np.meshgrid(points_x,points_y)

# points = np.tile(points,ny)
# points = np.reshape(points,(ny,nx))

# ------------------------ Show mesh -------------------------------------------
lim = max([w,L])
if(show_plots):
    plt.figure()
    plt.plot([0,L,L,0,0],[0,0,w,w,0])       #gives outline of the bar (what we are simulating our 1D conduction)
    plt.scatter(points_x,points_y,1,color='red')
    plt.ylim([-lim/10-lim,lim+lim/10])
    plt.xlim([-lim/10, lim+lim/10])



# define our "a" and "b" array/matrix
ap = np.zeros_like(np.squeeze(points_x))
aW = np.zeros_like(ap)
aE = np.zeros_like(ap)
# aS = np.zeros_like(ap)
# aN = np.zeros_like(ap)
b = np.transpose(np.zeros(np.shape(ap)[0]))

S = 0    # heat generation per m
nk1 = 0
nk2 = 0
k = k1

for i in range(nx):

    if((i == 0 or i == nx-1) and method == "a"):
        delta_L = dL/2
    else:
        delta_L = dL
        # A_c could be caluclated here as a function of the rod diameter.

    # set aW. Left most point is 0.
    if(i > 0):
        aW[i] = k*A_c/dw
    else:
        aW[i] = 0
    # set aE. Right most point is 0.
    if(i < nx-1):
        aE[i] = k*A_c/dw
    else:
        aE[i] = 0
    # if(j < ny-1):
    #     if(method == "a"):
    #         if(j > 0):
    #             aN[j,i] = k*dw/dL
    #         else:
    #             aN[j,i] = k*dw/(2*dL)
    #     else:
    #         aN[j,i] = k*dL/dw
    # else:
    #     aN[j,i] = 0
    # if(j > 0):
    #     if(method == "a"):
    #         if(j < ny-1):
    #             aS[j,i] = k*dw/dL
    #         else:
    #             aS[j,i] = k*dw/(2*dL)
    #     else:
    #         aS[j,i] = k*dw/dL
    # else:
    #     aS[j,i] = 0

    ap[i] = S/(delta_L) + (aW[i] + aE[i])      # ap - aw - ae - S = 0
        # ap[j,i] = S/(delta_L*delta_w) + (aW[j,i] + aE[j,i] + aS[j,i] + aN[j,i]    # 2D?

# make "A" matrix  (A*T = b)
A = np.zeros([nx,nx])

# remember, ap = aw + ae + an + as + S, so our aw and ae have to negative in the A matrix.
for i in range(nx):
    if(i > 0):
        A[i-1,i] = -aW[i]   
    if(i < nx-1):
        A[i+1,i] = -aE[i]
    A[i,i] = ap[i]

# set boundary conditions in A matrix and b vector
# make "b" vector, set with boundary conditions. (A*T = b)
b = np.zeros([nx,1])

# Dirichlet boundary conditions -> prescibed temperature at one end
if(b1_method == "temp"):
    A[:,0] = 0
    A[0,0] = 1          # Dirichlet boundary condition at x=0
    b[0] = T0           # Dirichlet boundary condition at x=0
if(b2_method == "temp"):
    A[:,-1] = 0
    A[-1,-1] = 1        # Dirichlet boundary condition at x=L
    b[-1] = TL          # Dirichlet boundary condition at x=L
# Neumann boundary conditions -> prescibed heat flux or insulated (q = 0)
if(b1_method == "insulated"):
    A[0,0] = -k/(dL/2)
    A[0,1]= k/(dL/2)
    b[0] = 0            # insulated = 0, prescribed heat flux, b = q0
elif(b1_method == "flux"):
    A[0,0] = -k/(dL/2)
    A[0,1]= k/(dL/2)
    b[0] = -qL

if(b2_method == "insulated"):
    A[-1,-1] = -k/(dL/2)
    A[-1,-2]= k/(dL/2)
    b[-1] = 0            # insulated = 0, prescribed heat flux, b = q0
elif(b2_method == "flux"):
    A[-1,-1] = -k/(dL/2)
    A[-1,-2]= k/(dL/2)
    b[-1] = -qR

# Robin boundary conditions


A = np.transpose(A)
T = la.solve(A,b)

print(T)

if(show_plots == True):
    plt.figure()
    try:
        plt.plot(np.transpose(points_x),T)
    except:
        plt.plot(points_x,T)
    plt.xlabel('Location (m)')
    plt.ylabel('Temp ($^\circ$C)')
    plt.show()

    contour_map_y = np.linspace(w/2,-w/2,ny)

    # define x,y coordinates
    plot_plane_x = points_x
    plot_plane_y = contour_map_y

    Y,X = np.meshgrid(plot_plane_y,plot_plane_x)
    reps = np.shape(plot_plane_y)
    T_plot = np.tile(np.array([T]),reps)
    # ny = np.size(contour_map_y)
    # nx = n
    T_plot = np.reshape(T_plot, (ny,nx))

    plt.figure()
    if(nx > 1 and ny > 1):
        plt.contourf(X,Y,T_plot)
    else:
        # T_plot_2D = T_plot.reshape(1,-1)
        plt.imshow(T_plot, aspect='auto', cmap='viridis',extent=[points_x.min(), points_x.max(), -w/2, w/2])
        plt.colorbar(label='Temp ($^\circ$C)')
    plt.ylim([-lim/10-lim,lim+lim/10])
    plt.xlim([0-dL,L+dL])
    plt.show()