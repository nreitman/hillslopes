#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Python2 code for exercise 4: hillslope evolution
written by Nadine on 2/8/2018

V2: REDO WITH SIMPLE ARRAYS TO DEBUG
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import os

degree = u"\u00b0" # unicode symbol for degree for labeling plots nicely

#################################################################################
#%% PART 2 - PLANAR HILLSLOPE WITH STEADILY UPLIFT VERTICAL FAULT
#################################################################################
# define variables
#kappa = .06
rhos = 1.33 * 1e3 #1.33 * 1e6 # density soil 1.33 g/cm3 --> [g/m^3]
rhor = 2.65 * 1e3 # 2.65 * 1e6 # density rock 2.65 g/cm3  --> [g/m^3]
#k = kappa * rhos
D = 50 * 1e-4 # [m^2/yr] from Fernandes and Dietrich 1997 in WRR
k = D/rhos

L = 300. # length of slope [meters]
dx = .5 # x step [meters]
x = np.arange(0,L+dx,dx) # initialize x array [meters]

dt = 1. # timestep [years]
tmax = 11000. # length of time for model run [years]
time = np.arange(0,tmax+dt+dt,dt) # initialize time array [years]
ideal_dt = (dx**2)/(2*D)
print('ideal timestep is less than: '+str(ideal_dt))
print('timestep is currently: '+str(dt)+' years')

H = np.ndarray(shape=(len(x),len(time)),dtype=float) # initialize array for soil thickness [meters]
H_initial = 1.5 # initial soil thickness - 1 meter
H[:,0] = H_initial   # start model at inital H [meters]
Hstar = .5 # scaling parameter for weathering rate. No idea what value to use. Starting with 1. 

Wo = 1e-6 #1e-5 # initial weathering rate mm/yr [meters/yr]
W = np.empty(shape=(len(x),len(time)),dtype=float) # initialize array to hold weathering rates
W[:,0] = Wo

z = np.zeros(shape=(len(x),len(time)),dtype=float) # intialize array for topo elevation [z(x,t)] [meters], add two b/c one cell on each end is for river height boundary condition
ztop = 175.
z_initial = np.linspace(ztop,0,num=len(x)) # initial topography [meters] - planar hillslope
z[:,0] = z_initial     # topo elevation for first timestep

fault = 150
line = z[np.int(fault/dx),0]/ztop
dip = 90.
alpha = np.int(np.degrees(np.arctan(ztop/L)))
slip_vertical = 1 # all slip is vertical
slip_horizontal = 0 # no slip is horizontal

Useismic = 3. # rate of fault uplift 1 m [m/yr] during seismic event
Uinterseis = 0. # rate of fault uplift during interseismic periods
U = np.zeros(len(time))
#U[10] = Useismic
for i in range(1,len(time)):
    if i % 2000 == 0:
        U[i] = Useismic
    elif i % 2000 != 0:
        U[i] = Uinterseis

plt.plot(time,U)
plt.xlabel('time [years]')
plt.ylabel('uplift rate [m/yr]')
plt.title('Uplift Through Time')
plt.show()


zb = np.zeros(shape=(len(x),len(time)),dtype=float) # initialize array for bedrock elevation
zb[:,0] = z[:,0] - H[:,0] # bedrock elevation [meters] - set to z_initial - H (soil thickness) = 1 meter
#zb[:,0] = z[:,0]

# initialize the rest of the arrays
Q = np.zeros(shape=(len(x)-1,len(time)),dtype=float)
dzdx = np.zeros(shape=(len(x)-1,len(time)),dtype=float)
dQdx = np.zeros(shape=(len(x)-2,len(time)),dtype=float)
dHdt = np.zeros(shape=(len(x)-2,len(time)),dtype=float)

# plot intial z (topo elev) and zb (bedrock elev)
plt.figure(figsize=(6,4))
plt.plot(x,zb[:,0],'brown',linestyle='--')
plt.plot(x,z[:,0],'mediumseagreen',linestyle='--')
plt.axvline(150,0,line,color='k',linewidth=0.85) # plot line at X = 150
plt.grid(color='lightgray',linestyle='--')
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('finite diff solution - initial conditions')
plt.text(255,1,str(alpha)+degree,color='blue')
plt.text(fault+1,z[np.int(fault/dx),0]+2,'dip = '+str(dip)+degree,color='black')
#plt.ylim(0,25)
plt.show()


#%% run finite diff loop
#plt.figure(figsize=(6,4))
#plt.plot(x,z_ss,'mediumseagreen',linestyle='--')
#plt.ylim(0,25)
#plt.grid(color='lightgray',linestyle='--')
#plt.xlabel('distance [m]')
#plt.ylabel('elevation [m]')
#plt.title('finite diff - steady state hillslope')

for i in range(0,len(time)-1):
    W[:,i] = Wo * (np.exp(-H[:,i]/Hstar))              # calculate weathering rate at this time i 
    dzdx[:,i] = np.diff(z[:,i]) / dx                    # calculate slope gradient at this time i
    Q[:,i] = - D * dzdx[:,i]                            # calculate soil flux for all x at this time i
    dQdx[:,i] = np.diff(Q[:,i]) / dx                    # calculate gradient in Q (soil flux) for all x at this time i
    dHdt[:,i] = (W[1:-1,i] * (rhor/rhos)) - dQdx[:,i] 
    dHdt[:,i] = (Wo * (rhor/rhos)) - dQdx[:,i]# calculate rate of change of soil thickness for all x at this time i 
    H[1:-1,i+1] = H[1:-1,i] + dHdt[:,i] * dt                        # update soil thickness at all x at next time i+1
    zb[0:np.int(fault/dx),i+1] = zb[0:np.int(fault/dx),i] + (U[i] * dt) # update bedrock elevation for upper plate at next time i+1 at tectonic rate
    zb[np.int(fault/dx):-1,i+1] = zb[np.int(fault/dx):-1,i]             # update bedrock elevation for lower plate at next time i+1 (doesn't change tectonically)
    z[:,i+1] = zb[:,i+1]  + H[:,i+1]                                              # update topo elevation for all x at next time i+1
    
    #if i % plots == 0:                          # plot at every plots timestep
       #plt.plot(x,z[:,i])                       # plot z (topo) for all x at time i
       #plt.text(-20, 12.5, 'time (years):'+ str((i))) # add label total time
       #plt.savefig('tmp'+str(i/plots)+'.png',bbox_inches="tight",dpi=150) # save plot for movie
 
#plt.text(-49, 23, 'time (years): '+ str((i))) # add label total time
#plt.show()


#%% plot output

plots = 200
       
plt.figure(figsize=(6,4))
#plt.plot(x,z_ss,'mediumseagreen',linestyle='--')
plt.ylim(85,105)
plt.xlim(140,160)
#plt.xlim(140,160)
plt.grid(color='lightgray',linestyle='--')
plt.axvline(150,0,line,color='k',linewidth=0.85) # plot line at X = 150
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('finite diff - tectonically uplifting hillslope')


for i in range(0,len(time)-1):#0,1,9,100,200,300:
    if i % plots == 0:
#range(0,len(time)-1):
#1,10,100,1000, 2000, 3000, 4000, 5000:
# range(0,int((tmax/dt)/100)):
#0,1,10000, 20000, 30000, 40000, 50000:
#1,10000, 20000, 30000, 40000, 50000:
#1,1e6, 2e6, 3e6, 4e6, 5e6:
#range(np.linspace(0,int(tmax/dt),num=6,dtype=int)):
        plt.plot(x,z[:,i],label=str(i*dt)+' years') 
        #plt.text(-49, 9.5, 'time (years): '+ str((i*dt))) # add label total time
        #plt.legend()
plt.text(155,115, 'time (years): '+ str(int((i*dt)))) # add label total time
#plt.text(150,31,'Wo: '+str(Wo)+' mm/yr')
#plt.text(0,8,'Uplift: '+str(Useismic)+' mm/yr')