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

#################################################################################
#%% ANALYTICAL SOLUTION - STEADY STATE CONVEX UP HILLSLOPE WITH INCISING STREAM
#################################################################################
# define variables
U = 1 * .001 # rate of baselevel fall 1 mm.yr [meters/yr]
L = 50 # length of slope [meters]
#k =  10000 # ????
rhos = 1.33 * 1e6 # density soil 1.33 g/cm3 --> [g/m^3]
rhor = 2.65 * 1e6 # density rock 2.65 g/cm3  --> [g/m^3]
#kappa = k / rhos
kappa = .06 # [meters^2/yr] (value from geomechanics notes 10/8)
dx = 1. # x step [meters]
x = np.arange(-L,L+dx,dx) # initialize x array

# calculate topography (z)
z = (U / (2 * kappa)) * ((L**2) - (x**2))

# plot 
plt.figure(figsize=(6,4))
plt.plot(x,z,'mediumseagreen')
plt.ylim(0,25)
plt.grid(color='lightgray',linestyle='--')
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('analytical solution - steady state hillslope')
plt.show()

z_ss = z # save steady state z profile as initial condition z topo for finite diff solution

#################################################################################
#%% PART 1 - CONVEX UP HILLSLOPE WITH INCISING STREAM
#################################################################################
# define variables
plots = 1

#kappa = .06
rhos = 1.33 * 1e6 # density soil 1.33 g/cm3 --> [g/m^3]
rhor = 2.65 * 1e6 # density rock 2.65 g/cm3  --> [g/m^3]
#k = kappa * rhos
D = 50 * 1e-4 # [m^2/yr] from Fernandes and Dietrich 1997 in WRR
k = D/rhos

L = 50. # length of slope [meters]
dx = .5 # x step [meters]
x = np.arange(-L,L+dx,dx) # initialize x array [meters]

dt = 10. # timestep [years]
tmax = 5000000. # length of time for model run [years]
time = np.arange(0,tmax+dt+dt,dt) # initialize time array [years]
ideal_dt = (dx**2)/(2*D)
print('ideal timestep is less than: '+str(ideal_dt))
print('timestep is currently: '+str(dt)+' years')

H = np.ndarray(shape=(len(x),len(time)),dtype=float) # initialize array for soil thickness [meters]
H_initial = 1. # initial soil thickness - 1 meter
H[:,0] = H_initial   # start model at inital H [meters]
#H[0,:] = 0. # soil thickness 0 at stream channels for all of time
#H[-1,:] = 0. # soil thickness 0 at stream channels for all of time
Hstar = .5 # scaling parameter for weathering rate. No idea what value to use. Starting with 1. 

Wo = 1e-3 #1e-5 # initial weathering rate mm/yr [meters/yr]
W = np.empty(shape=(len(x),len(time)),dtype=float) # initialize array to hold weathering rates
W[:,0] = Wo

U = 1e-5 # rate of baselevel fall (or uplift rate) 1 mm/yr [meters/yr] doesn't change in time or space...YET

z = np.zeros(shape=(len(x),len(time)),dtype=float) # intialize array for topo elevation [z(x,t)] [meters], add two b/c one cell on each end is for river height boundary condition
z_initial = z_ss # initial topography [meters] - output from steady state solution
z[:,0] = 1.     # topo elevation for first timestep, except at channels, set below
z[0,:] = 0.     # channel height = 0 for all time on left side (boundary condition)
z[-1,:] = 0.    # channel height = 0 for all time on right side (boundary condition)

zb = np.zeros(shape=(len(x),len(time)),dtype=float) # initialize array for bedrock elevation
zb[:,0] = z[:,0] - H[:,0] # bedrock elevation [meters] - set to steady state solution - H (soil thickness) = 1 meter
zb[0,:] = 0.    # bedrock elevation = 0 for all of time at channels
zb[-1,:] = 0.   # bedrock elevation = 0 for all of time at channels

# initialize the rest of the arrays
Q = np.zeros(shape=(len(x)-1,len(time)),dtype=float)
dzdx = np.zeros(shape=(len(x)-1,len(time)),dtype=float)
dQdx = np.zeros(shape=(len(x)-2,len(time)),dtype=float)
dHdt = np.zeros(shape=(len(x)-2,len(time)),dtype=float)

# plot intial z (topo elev) and zb (bedrock elev)
plt.figure(figsize=(6,4))
plt.plot(x,zb[:,0],'brown',linestyle='--')
plt.plot(x,z[:,0],'mediumseagreen',linestyle='--')
#plt.ylim(0,25)
plt.grid(color='lightgray',linestyle='--')
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('finite diff solution - initial conditions')
plt.show()


#%% run finite diff loop
#plt.figure(figsize=(6,4))
#plt.plot(x,z_ss,'mediumseagreen',linestyle='--')
#plt.ylim(0,25)
#plt.grid(color='lightgray',linestyle='--')
#plt.xlabel('distance [m]')
#plt.ylabel('elevation [m]')
#plt.title('finite diff - steady state hillslope')

for i in range(len(time)-1):
    W[:,i] = Wo * (np.exp(-H[:,i]/Hstar))               # calculate weathering rate at this time i 
    dzdx[:,i] = np.diff(z[:,i]) / dx                    # calculate slope gradient at this time i
    Q[:,i] = - D * dzdx[:,i]                            # calculate soil flux for all x at this time i
    dQdx[:,i] = np.diff(Q[:,i]) / dx                    # calculate gradient in Q (soil flux) for all x at this time i
    dHdt[:,i] = (W[1:-1,i] * (rhor/rhos)) - dQdx[:,i]   # calculate rate of change of soil thickness for all x at this time i
    H[1:-1,i+1] = H[1:-1,i] + dHdt[:,i] * dt                        # update soil thickness at all x at next time i+1
    zb[1:-1,i+1] = zb[1:-1,i] + (-W[1:-1,i] * dt) + (U * dt)          # update bedrock elevation for all x at next time i+1
    z[1:-1,i+1] = zb[1:-1,i+1] + H[1:-1,i+1]            # update topo elevation for all x at next time i+1
    
    #if i % plots == 0:                          # plot at every plots timestep
       #plt.plot(x,z[:,i])                       # plot z (topo) for all x at time i
       #plt.text(-20, 12.5, 'time (years):'+ str((i))) # add label total time
       #plt.savefig('tmp'+str(i/plots)+'.png',bbox_inches="tight",dpi=150) # save plot for movie
 
#plt.text(-49, 23, 'time (years): '+ str((i))) # add label total time
#plt.show()


#%% plot output
plt.figure(figsize=(6,4))
#plt.plot(x,z_ss,'mediumseagreen',linestyle='--')
#plt.ylim(15,25)
plt.ylim(0,10)
plt.grid(color='lightgray',linestyle='--')
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('finite diff - steady state hillslope')

for i in 1,10,100,1000, 2000, 3000, 4000, 5000:
# range(0,int((tmax/dt)/100)):
#0,1,10000, 20000, 30000, 40000, 50000:
#1,10000, 20000, 30000, 40000, 50000:
#1,1e6, 2e6, 3e6, 4e6, 5e6:
#range(np.linspace(0,int(tmax/dt),num=6,dtype=int)):
    plt.plot(x,z[:,i],label=str(i*dt)+' years') 
    plt.legend()
