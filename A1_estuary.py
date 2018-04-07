# coding: utf-8
#Do Balance non-dim
import numpy as np
import math
import pylab as py

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

def estuary(scal, width):
    #set up arrays (python)
    x=np.zeros((nspace),dtype=float) #space domain
    C=np.ones((mtime,nspace),dtype=float) #DO
    L=np.zeros((mtime,nspace),dtype=float) #CBOD
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
    print(delt_sec*mtime/3600/24) 
    #Streeter Phelps
    Csp = np.zeros((nspace),dtype=float) #space domain
    fac = Kr*Lin/(Ka-Kr)
    for ii in range(0,nspace):
        x[ii]=ii*delx #defining range x
        #Csp[ii] = Cs-fac*(math.exp(-Kr*x[ii])-math.exp(-Ka*x[ii]))
        Csp[ii]=1-fac*(math.exp(-Kr*x[ii])-math.exp(-Ka*x[ii]))
    
    # BOD at outfall
    L[0,0]=Lin
    for tt in range(0,mtime-1):
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
    xscale=(D/u)/1000 #km
    return x, C, Csp, xscale

x, C, Csp, xscale = estuary(.1, 5)
#x2, C2, Csp2, xscale2 = estuary(.1, 50)
#x3, C3, Csp3, xscale3 = estuary(0.1, 200)
#x4, C4, Csp4, xscale4 = estuary(0.3, width)
#x5, C5, Csp5, xscale5 = estuary(0.3, 50)
#x6, C6, Csp6, xscale6 = estuary(0.001, 50)
py.figure(1)

py.plot(xscale*x,C[mtime-1,:],'.r', label='Run 1')
py.plot(xscale*x,Csp,'r', label='SS Streeter-Phelps')
#py.plot(xscale2*x2,C2[mtime-1,:],'--y', label='Run 2')
#py.plot(xscale2*x2,Csp2,'y', label='SS Streeter-Phelps')
#py.plot(xscale3*x3,C3[mtime-1,:],'.y', label='Run 3')
#py.plot(xscale3*x3,Csp3,'y', label='SS Streeter-Phelps')
#py.plot(xscale4*x4,C4[mtime-1,:],'.b', label='Run 4')
#py.plot(xscale4*x4,Csp4,'b', label='SS Streeter-Phelps')
#py.plot(xscale5*x5,C5[mtime-1,:],'.g', label='Run 5')
#py.plot(xscale5*x5,Csp5,'g', label='SS Streeter-Phelps')
#py.plot(xscale6*x6,C6[mtime-1,:],'--r', label='Run 6')
#py.plot(xscale6*x6,Csp6,'r', label='SS Streeter-Phelps')
py.xlabel('distance downstream (m)')
py.ylabel('concentration DO (kg/m^3)')
py.title('Disolved Oxygen being dispersed in a tidal estuary')
py.legend()
py.show()
