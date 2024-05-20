import math
import numpy as np
import sys
import glob
import re
import os
from subprocess import call
import matplotlib.pyplot as plt



#Plot settings
label_size = 30
tick_size = 18



#Loading grid and time data
inpt = np.loadtxt("Output/exit_charge_species.dat")


#Removing duplicates and ordering
u, ind = np.unique(inpt[:,0], return_index=True)

nind = np.shape(ind)[0]

data = inpt[ind,:]


#Other inputs
T = np.loadtxt("Output/time.dat")
X = np.loadtxt("Output/xgrid.dat")
Y = np.loadtxt("Output/ygrid.dat")


nt = T[0]
pint = T[1]
dt = T[2]



xwidth = X[-1] - X[0]


#Getting time information
time = data[:,0]*dt
dpt = time[1] - time[0]

#Computing Nec1
dN_ec1 = data[:,3] - data[:,4] - data[:,1]


#Clumping parameter
CQ = 5.0/3.0*1.0e6


#Converting to fluxes
gamma_ec1 = CQ*dN_ec1/dpt/xwidth

#Obtaining jec1
j_ec1 = gamma_ec1*1.6e-19



#Plotting the exiting number of particles against time - normalized against jM
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

plt.plot(time*1.0e6,j_ec1/4.0,'g',linewidth=2)

plt.xlabel(r"$t$ $(\mu s)$", fontsize=label_size)
plt.ylabel(r"$j_{ec1} / j_M$ $\%$", fontsize=label_size)

#plt.axis('equal')

plt.xlim([0.0,time[-1]*1.0e6])
plt.ylim([0.0,200.0])


#Shifting axis slightly
pos1 = ax.get_position()
pos2 = [pos1.x0+0.045,pos1.y0+0.02,pos1.width,pos1.height]
ax.set_position(pos2)

plt.tick_params(labelsize=tick_size)

plt.savefig("exit_current.png")

plt.close(fig)




#Saving the data
OUT = np.zeros([np.shape(time)[0],2])

OUT[:,0] = time*1.0e6
OUT[:,1] = j_ec1/4.0

np.savetxt("jec1_us_normalized.dat",OUT,fmt='%.8e')

















