#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Python2 code for exercise 4: hillslope evolution
written by Nadine on 2/8/2018
part 1 = CONVEX UP HILLSLOPE WITH INCISING STREAM
part 2 = PLANAR HILLSLOPE WITH INTERMITTENT FAULT
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
k =  10000 # ????
rhos = 1.33 * 1e6 # density soil 1.33 g/cm3 --> [g/m^3]
rhor = 2.65 * 1e6 # density rock 2.65 g/cm3  --> [g/m^3]
kappa = k / rhos
kappa = .06 # [meters^2/yr] (value from geomechanics notes 10/8)
dx = 1 # x step [meters]
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

#################################################################################
#%% PART 1 - CONVEX UP HILLSLOPE WITH INCISING STREAM
#################################################################################
# define variables
plots = 100

kappa = .06
rhos = 1.33 * 1e6 # density soil 1.33 g/cm3 --> [g/m^3]
rhor = 2.65 * 1e6 # density rock 2.65 g/cm3  --> [g/m^3]
k = kappa/rhos

L = 50. # length of slope [meters]
dx = 1. # x step [meters]
x = np.arange(-L,L+dx,dx) # initialize x array [meters]

dt = 10. # timestep [years]
tmax = 1000000. # length of time for model run [years]
time = np.arange(0,tmax+dt,dt) # initialize time array [years]

H = np.ndarray(shape=(len(x),len(time)),dtype=float) # initialize array for soil thickness [meters]
H_initial = 1 # initial soil thickness - 1 meter
H[:,0] = H_initial   # start model at inital H [meters]
H[0,:] = 0 # soild thickness 0 at stream channels for all of time
H[-1,:] = 0 # soild thickness 0 at stream channels for all of time
Hstar = 1 # scaling parameter for weathering rate. No idea what value to use. Starting with 1. 

Wo = 2 * .001 # initial weathering rate 2 mm/yr [meters/yr]
W = np.zeros(len(time)) # initialize array to hold weathering at each timestep

U = 1 * .001 # rate of baselevel fall (or uplift rate) 1 mm.yr [meters/yr] doesn't change in time or space...YET

z = np.ndarray(shape=(len(x),len(time)),dtype=float) # intialize array for topo elevation [z(x,t)] [meters], add two b/c one cell on each end is for river height boundary condition
z_initial = 1 # initial topography [meters]
z[:,0] = z_initial # topo elevation for first timestep is 1 meter, except at channels, set below
z[0,:] = 0  # channel height = 0 for all time on left side (boundary condition)
z[-1,:] = 0 # channel height = 0 for all time on right side (boundary condition)

zb = np.ndarray(shape=(len(x),len(time)),dtype=float) # initialize array for bedrock elevation
zb[:,0] = 0 # bedrock elevation [meters] - set to zero for the first timestep
zb[0,:] = 0 # bedrock elevation = 0 for all of time at channels
zb[-1,:] = 0 # bedrock elevation = 0 for all of time at channels

W = np.ndarray(shape=(len(x),len(time)),dtype=float)
Q = np.ndarray(shape=(len(x)-1,len(time)),dtype=float)
dzdx = np.ndarray(shape=(len(x)-1,len(time)),dtype=float)
dQdx = np.ndarray(shape=(len(x)-2,len(time)),dtype=float)
dHdt = np.ndarray(shape=(len(x)-2,len(time)),dtype=float)
#%%
plt.figure(figsize=(6,4))
plt.plot(x,z,'mediumseagreen')
plt.ylim(0,25)
plt.grid(color='lightgray',linestyle='--')
plt.xlabel('distance [m]')
plt.ylabel('elevation [m]')
plt.title('finite diff - steady state hillslope')

for i in range(len(time)-1):
    W[:,i] = Wo * (np.exp(-H[:,i]/Hstar))       # calculate weathering rate for time i. 
    dzdx[:,i] = np.diff(z[:,i]) / dx            # calculate slope gradient
    Q[:,i] = - k * dzdx[:,i]                    # calculate soil flux for all x at time i
    dQdx[:,i] = np.diff(Q[:,i]) / dx            # calculate gradient in Q (soil flux) for all x at time i
    dHdt[:,i] = (W[1:-1,i] * (rhor/rhos)) - dQdx[:,i] # calculate rate of change of soil thickness for all x at time i
    H[1:-1,i] = dHdt[:,i] * dt                     # update soil thickness at all x for time i
    zb[1:-1,i] = (W[1:-1,i] * dt) + (U * dt)            # update bedrock elevation for all x at time i
    z[1:-1,i+1] = zb[1:-1,i] + H[1:-1,i]                # update topo elevation for all x at this time i
    
    if i % plots == 0:                          # plot at every plots timestep
       plt.plot(x,z[:,i])                       # plot z (topo) for all x at time i
       #plt.text(-20, 12.5, 'time (years):'+ str((i))) # add label total time
       #plt.savefig('tmp'+str(i/plots)+'.png',bbox_inches="tight",dpi=150) # save plot for movie
 
plt.text(-20, 12.5, 'time (years): '+ str((i))) # add label total time
plt.show()


