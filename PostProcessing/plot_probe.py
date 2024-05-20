import math
import numpy as np
import sys
import glob
import re
import os
from subprocess import call
import matplotlib.pyplot as plt
import warnings

#Setting the working directory to the top level code folder
script_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = script_dir.split('/')

if script_dir[-1] == "PostProcessing":
	script_dir = script_dir[:-1]

work_dir = '/'.join(script_dir)
os.chdir(work_dir)


### START USER INPUT ###

#Directory in Output/ to search for data
data_dir = "ProbeDiagnostics"

#Number of the phase space diagnostic corresponding to the number provided in your LTP-PIC input file
diag_num = 1

#Directory to create for saving time dependent plots. Time averaged plots are saved in the top directory
print_dir = "Probe" + str(diag_num)

#Plot lower bound in t-direction. If set to 'None' it will be chosen automatically.
t_min = None

#Plot upper bound in t-direction. If set to 'None' it will be chosen automatically.
t_max = None

#Plot lower bound in y-direction. If set to 'None' it will be chosen automatically.
y_min = None

#Plot upper bound in y-direction. If set to 'None' it will be chosen automatically.
y_max = None

### END USER INPUT ###

#First check if the probe exists
if (os.path.exists("Output/" + data_dir + "/p%d.dat" % (diag_num)) == False):
	print("Probe number %d does not exist\nEXITING\n")
	sys.exit()

#Getting the variable name
with open("Output/" + data_dir + "/p%d.dat" % (diag_num)) as f:
    var = f.readline().strip('\n')

#Getting the data
DAT = np.loadtxt("Output/" + data_dir + "/p%d.dat" % (diag_num), skiprows=1)

T = DAT[:,0]
DAT = DAT[:,1]

#print(var)
#print(T)
#print(DAT)


#Check if the probe print directory exists. Otherwise create it.
print_dir = "ProbePlots"

if (os.path.isdir(print_dir) == False):
	call(['mkdir',print_dir])



#Plot settings
label_size = 30
tick_size = 18

#Shifting axis
xshift = 0.06
yshift = 0.02

#Turning off numpy warnings if an empy data set is loaded
warnings.simplefilter("ignore")



#Setting the y-variable name
if var == 'density':
	ylab = r"$Density$ $(1/m^3)$"
elif var == 'chargedensity':
	ylab = r"$Charge Density$ $(C/m^3)$"
elif var == 'potential':
	ylab = r"$Potential$ $(V)$"
elif var == 'electricfieldx':
	ylab = r"$E_x$ $(V/m)$"
elif var == 'electricfieldy':
	ylab = r"$E_y$ $(V/m)$"
elif var == 'electricfieldz':
	ylab = r"$E_z$ $(V/m)$"
else:
	ylab = var



#Setting to user defined min/max if defined
if (t_min == None):
	tmin = T[0]
else:
	tmin = t_min

if (t_max == None):
	tmax = T[-1]
else:
	tmax = t_max

if (y_min == None):
	ymin = np.min(DAT)
else:
	ymin = y_min

if (y_max == None):
	ymax = np.max(DAT)
else:
	ymax = y_max



#Adjusting the limits for velocity to give some buffer space
	yh = ymax - ymin
	ymax += 0.2*yh
	#ymin -= 0.2*yh


#Plotting
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

plt.plot(T, DAT, linewidth=2)

plt.xlabel(r"$t$ $(s)$", fontsize=label_size)
plt.ylabel(ylab, fontsize=label_size)

plt.xlim([tmin,tmax])
plt.ylim([ymin,ymax])

plt.tick_params(labelsize=tick_size)


#Shifting axis slightly
pos1 = ax.get_position()
pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
ax.set_position(pos2)

plt.savefig(print_dir + "/Probe%d.png" % diag_num)

plt.close(fig)
