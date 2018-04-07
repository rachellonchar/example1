#8561 Assignment 1
#advection dispersion 
#point source

import numpy as np
import math
import pylab as plt
import sys

def C_reaction(D, k):
    #set up space-time grid
    nspace = 100 #space steps 0,1,2.....,98,nspace-1
    mtime = 100 #time steps 0,1,2.....,98,mtime-1 
    #fixed parameters
    u = 1. #velocity
    if k == 0:
        kcon = 0
    else:
        kcon = 0.001
    delx = 1*D*u #ensures Courant condition
    delt = .25*delx*delx/D #ensures Fourier condition
    #set up arrays (python)
    x = np.zeros((nspace),dtype=float)
    C = np.zeros((mtime,nspace),dtype=float)
    C_analytical = np.zeros((mtime,nspace),dtype=float)
    for ii in range(0,nspace): 
         x[ii] = ii*delx
    #boundary condition at x = 0 and time t = 0
    C[0,0] = 1
    C_analytical[0,0] = 1
    for tt in range(0,mtime-1):
        C[tt+1,0] = 1 #point source
        C_analytical[tt+1,0] = 1 #point source
        for ii in range(1,nspace-1): #go one short to keep C[tt,nsapce-1]=0
            #numerical solution:
            ap = 1 - (delt/delx/delx)*(D+D) - delt*k
            aw = (delt/delx/delx)*D + (u*delt/2./delx)
            ae = (delt/delx/delx)*D - (u*delt/2./delx)
            C[tt+1,ii] = ap*C[tt,ii] + aw*C[tt,ii-1] + ae*C[tt,ii+1] + delt*kcon
            if k == 0: #conservative
                #analytical solution:
                erfc_m = math.erfc((delx*ii - u*delt*(tt+1))/(2*(D*delt*(tt+1))**0.5))
                erfc_p = math.erfc((delx*ii + u*delt*(tt+1))/(2*(D*delt*(tt+1))**0.5))
                C_analytical[tt+1, ii] = 1/2*(erfc_m + np.exp(u*delx*ii/D)*erfc_p)
            else:
                #analytical upper bound:
                C_analytical[tt+1, ii] = np.exp(-k*ii*delx/u)
    return x, C, C_analytical

#parameters to plot
D = [1, 2] #dispersion
k = [0.4, 0.01] #reaction rate, first order
mtime = 100

#1a.
#conservative tracer
fig1a = plt.figure()
x, C1, C1_analytical = C_reaction(D[0], 0)
x, C2, C2_analytical = C_reaction(D[1], 0)
plt.figure(1)
plt.plot(x,C1[mtime-1,:], 'bo', label='numerical solution (D=1)')
plt.plot(x,C1_analytical[mtime-1,:], 'b', label='analytical solution (D=1)')
plt.plot(x,C2[mtime-1,:], 'ro', label='numerical solution (D=2)')
plt.plot(x,C2_analytical[mtime-1,:], 'r', label='analytical solution (D=2)')
plt.xlabel('distance downstream')
plt.ylabel('concentration')
plt.title('Conservative tracer from a point source (x=0)')
plt.legend()
plt.savefig('g/hw1_1a.png')
plt.close(fig1a)

#1b.
#first order reaction
fig1b = plt.figure()
x, C14, C14_an = C_reaction(D[0], k[0])
x, C24, C24_an = C_reaction(D[1], k[0])
x, C101, C101_an = C_reaction(D[0], k[1])
x, C201, C201_an = C_reaction(D[1], k[1])
plt.figure(2)
plt.plot(x,C14[mtime-1,:], 'b', label='D = 1, k = 0.4')
plt.plot(x,C14_an[mtime-1,:], 'b--', label='upper bound')
plt.plot(x,C24[mtime-1,:], 'g', label='D = 2, k = 0.4')
plt.plot(x,C24_an[mtime-1,:], 'g--', label='upper bound')
plt.plot(x,C101[mtime-1,:], 'r', label='D = 1, k = 0.01')
plt.plot(x,C101_an[mtime-1,:], 'r--', label='upper bound')
plt.plot(x,C201[mtime-1,:], 'y', label='D = 2, k = 0.01')
plt.plot(x,C201_an[mtime-1,:], 'y--', label='upper bound')
plt.xlabel('distance downstream')
plt.ylabel('concentration')
plt.title('Point source tracer (x=0) with first order reaction')
plt.legend()
plt.savefig('g/hw1_1b.png')
plt.close(fig1b)
