# coding: utf-8
#Do Balance non-dim
import numpy as np
import math
import pylab as py
import pandas as pd
import sys

#fig = plt.figure() 

nspace=1000 #space steps 0,1,2.....,98,nspace-1
mtime=2000 #tiome steps 0,1,2.....,98,mtime-1

#Fixed items UNITS m, second, Kg, m^3

#stream conditions
u=1.0 # velocity m/s
#width=5.0 #m
h=0.5
Slope=0.001
ustar=math.sqrt(Slope*h*9.81) #shear velocity

#net photo-res Kg/m^\[LongDash]given
pnet=5.0/1000.0/3600.0/24 #net of photo-res kg/m^3-s
SB=0

#given initial conditions in water and effluent
Cs=8.0/1000.0 #saturation kg/m^3
Lin=0.6 #input CBOD kg/m^3
Lin = .1/Cs
#Cs = Cs/Cs
#Lin = Lin/Cs

def estuary(scal, width, days=2, param=0):
    #Can enhance value for more rapid calculation\[LongDash]explained in class
    Fac=5
    # Longitudinal Dispersion\[LongDash]Fisher 1979 m^2/s 0.011*u^2*w^2/(hu^*)
    Dreal=Fac*0.011*u*u*width*width/h/ustar
    # aeration rate constant (1/s) Wikipedia on page on Streeter-Phelps
    Ka=2.148*math.pow(u,0.878)/math.pow(h,1.48)/3600/24
    #CBOD reaction rat\[LongDash]set here as fraction of Ka
    Kr=scal*Ka
    Dalt=0.005*u*u/Kr
    D=max(Dreal,Dalt) #if adv dominated will use Max value of D possible
    #advection dominated if Re=DKr/(u^2) < .005
    Rey = D*Kr/(u**2)
    
    #dimensonless K's
    Kr=D*Kr/u/u
    Ka=D*Ka/u/u
    #delx=1*D/u #space step
    #delt=0.25*delx*delx/D #time step  
    #set space and time step
    delx=min(1,500*u/D)
    #Must be less than 2 BUT for accuracy physical size<=500m
    delt=.5/(2.0/delx/delx+Ka+Kr) 
    delt_sec = delt*D/(u**2)
    min_days = delt_sec*mtime/3600/24
    steps_days = days*3600*24/delt_sec
    if steps_days < mtime:
        tsteps = math.ceil(mtime)
        elapsed_days = min_days
    else:
        tsteps = math.ceil(steps_days)
        elapsed_days = delt_sec*tsteps/3600/24
    xscale=(D/u)/1000 #km
    #if only want these paramters and not to run the entire simulation
    if param == 1:
        return Ka, Kr, D, Rey, xscale
    elif param == 2:
        return tsteps, elapsed_days, delt_sec/3600/24
    else:
        #set up arrays (python)
        x=np.zeros((nspace),dtype=float) #space domain
        C=np.ones((tsteps,nspace),dtype=float) #DO
        L=np.zeros((tsteps,nspace),dtype=float) #CBOD
        #Streeter Phelps
        Csp = np.zeros((nspace),dtype=float) #space domain
        fac = Kr*Lin/(Ka-Kr)
        for ii in range(0,nspace):
            x[ii]=ii*delx #defining range x
            #Csp[ii] = Cs-fac*(math.exp(-Kr*x[ii])-math.exp(-Ka*x[ii]))
            Csp[ii]=1-fac*(math.exp(-Kr*x[ii])-math.exp(-Ka*x[ii]))
        
        # BOD at outfall
        L[0,0]=Lin
        for tt in range(0,tsteps-1):
            #update CBOD
            for ii in range(1,nspace-1):
                #regular fluxes
                qwest=(L[tt,ii-1]+L[tt,ii])/2
                qwest=qwest+(L[tt,ii-1]-L[tt,ii])/delx
                qeast=-(L[tt,ii+1]+L[tt,ii])/2
                qeast=qeast+(L[tt,ii+1]-L[tt,ii])/delx
                #Balance
                EXTRA=Kr*L[tt,ii]
                #check for anoxic conditions
                if C[tt,ii]==0:
                    EXTRA=Ka
                L[tt+1,ii]=L[tt,ii]+(delt/delx)*(qwest+qeast)-delt*EXTRA
            #set BOD at outfall
            L[tt+1,0]=Lin #upstream 
            L[tt+1,nspace-1]=L[tt+1,nspace-2] #downstream
            
            #update Do
            for ii in range(1,nspace-1):
                #regular fluxes
                qwest=(C[tt,ii-1]+C[tt,ii])/2
                qwest=qwest+(C[tt,ii-1]-C[tt,ii])/delx
                qeast=-(C[tt,ii+1]+C[tt,ii])/2
                qeast=qeast+(C[tt,ii+1]-C[tt,ii])/delx
                EXTRA=-Kr*L[tt,ii]+Ka*(1.0-C[tt,ii])
                C[tt+1,ii]=C[tt,ii]+(delt/delx)*(qwest+qeast)+delt*(EXTRA)
                #check for anoixic condition
                if C[tt+1,ii]<0:
                    C[tt+1,ii]=0
                else:
                    None
            #set down stream condition
            C[tt+1,nspace-1]=C[tt+1,nspace-2]
        return x, C, Csp

runs = [1, 2, 3, 4, 5, 6]
widths = [5, 50, 200, 5, 50, 50]
scals = [.1, .1, .1, .3, .3, 0.4]
locs = ["g/run1.png", "g/run2.png", "g/run3.png", "g/run4.png", "g/run5.png", "g/run6.png"]
locs_e = ["g/erun1.png", "g/erun2.png", "g/erun3.png", "g/erun4.png", "g/erun5.png", "g/erun6.png"]

def plott(run, error=1):
    run = run - 1 #adjust for indexing
    width = widths[run]
    scal = scals[run]
    loc = locs[run]
    loc2 = locs_e[run]
    Ka, Kr, D, R, xscale = estuary(scal, width, 1, 1)
    #print some parameters for the run
    if R < 0.005:
        dom = 'advection dominated'
    elif R > 20:
        dom = 'dispersion dominated'
    else:
        dom = 'combination of advection and dispersion'
    # print relevent data
    print('Run' + str(run+1) + ': width=' + str(width) + ' m, Kr=' + str(scal) + 'Ka')
    print("   Reynold's Number: " + str(R) + ', so ' + str(dom))
    print('D, Ka, Kr = ')
    print(D, Ka, Kr)
    fig = py.figure()
    #runs for 1, 2, and 4 years
    #steps, actual days, time step size (in days)
    s, d, dt = estuary(scal, width, 1, 2)
    s2, d2, dt2 = estuary(scal, width, 2, 2)
    s4, d4, dt4 = estuary(scal, width, 4, 2)
    print('---')
    print('(for intended days elapsed = 1, 2, 4)')
    print('actual days elapsed')
    print(d, d2, d4)
    print('step sizes:')
    print(dt, dt2, dt4)
    print('steps taken:')
    print(s, s2, s4)
    #run simulation
    x, C, Csp = estuary(scal, width, 1)
    x2, C2, Csp2 = estuary(scal, width, 2)
    x4, C4, Csp4 = estuary(scal, width, 4)
    diff1, diff2, diff4 = abs(C-Csp), abs(C2-Csp2), abs(C2-Csp2)
    diff = [diff1, diff2, diff4]
    #graphing
    py.plot(xscale*x,Csp,'r', label='1 day')
    py.plot(xscale*x,C[s-1,:],'.r')
    py.plot(xscale*x2,Csp2,'y', label='2 days')
    py.plot(xscale*x2,C2[s2-1,:],'.y')
    py.plot(xscale*x4,Csp4,'b', label='4 days')
    py.plot(xscale*x4,C4[s4-1,:],'.b')
    #return diff1, diff2, diff4
    #graph labels
    py.xlabel('distance downstream (km)')
    py.ylabel('concentration DO (unitless, 0=no DO, 1=saturation)')
    #py.xlabel('distance downstream (m)')
    #py.ylabel('concentration DO (kg/m^3)')
    py.title('DO in a tidal estuary, Run ' + str(run+1) + 
    ': width = ' + str(width)+ ' m, Kr = ' + str(scal) + 'Ka' +
    '\n(solid lines are the Streeter-Phelps steady state solutions)')
    py.legend() 
    #showing plots together
    py.savefig(loc)
    py.close(fig)
    if error == 1:
        fig2 = py.figure()
        py.plot(xscale*x,C[s-1,:]-Csp,'--r', label='1 day')
        py.plot(xscale*x2,C2[s2-1,:]-Csp2,'--y', label='2 days')
        py.plot(xscale*x4,C4[s4-1,:]-Csp4,'--b', label='4 days')
        #graph labels
        py.xlabel('distance downstream (km)')
        py.ylabel('difference in concentration, C[i] - C(Streeper-Phelps))')
        py.title('Difference between numerical and analytical methods, \nRun ' + str(run+1) + 
        ': width = ' + str(width)+ ' m, Kr = ' + str(scal) + 'Ka')
        py.legend() 
        #showing plots together
        #py.show()
        py.savefig(loc2)
        py.close(fig2)
    else:
        None
    print(' ')

#plott(1)
#plott(2)
plott(3)
plott(4)
plott(5)
#plott(6)

