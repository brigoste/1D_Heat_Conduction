import numpy as np
from numpy import linalg as la

def main(T0,TL,points,k,Area_vec,dx,Su):
    indecies = np.linspace(1,len(points)-2,1)   # only deal with interior nodes
    
    A = np.zeros((len(points),len(points)))
    
    # solve for ap, aw, ae, as
    for i in indecies:
        if(i == 1):
            ae = Area_vec[i+1]*k[i+1]/(dx)
            aw = Area_vec[i]*k[i]/(dx/2)
        elif(i == len(points)-2):
            ae = Area_vec[i+1]*k[i+1]/(dx/2)
            aw = Area_vec[i]*k[i]/(dx)
        else:
            ae = Area_vec[i+1]*k[i+1]/(dx)
            aw = Area_vec[i]*k[i]/(dx)

        S = Su*dx*((Area_vec[i]+Area_vec[i+1])/2)  # volumetric heat generation term
        ap = ae + aw - S
        
        # save ae, aw, and ap terms in A
        A[i,i-1] = -aw
        A[i,i+1] = -ae
        A[i,i] = ap

    #set temperature boundary conditions
    A[0,:] = 0
    A[-1,:] = 0
    A[0,0] = 1
    A[-1,-1] = 1

    b = np.zeros(len(points))
    b[0] = T0
    b[-1] = TL

    # solve the system of equations
    T = la.solve(A,b)

    return T