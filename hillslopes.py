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

#%% ANALYTICAL SOLUTION - STEADY STATE CONVEX UP HILLSLOPE WITH INCISING STREAM

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

#%% PART 1 - CONVEX UP HILLSLOPE WITH INCISING STREAM



#%% PART 2 - PLANAR HILLSLOPE WITH INTERMITTENT FAULT

