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
data_dir = "PhaseDiagnostics"

#Number of the phase space diagnostic corresponding to the number provided in your LTP-PIC input file
diag_num = 1

#Directory to create for saving time dependent plots. Time averaged plots are saved in the top directory
print_dir = "PhaseDiagnostic" + str(diag_num) + "_plots"

#Start time (s)
tstart = 0.0

#End time (s)
tend = 1.0

#Select the variable to plot on the x-axis
variable_x = 'x'

#Select the variable to plot on the x-axis
variable_y = 'y'

#Plot lower bound in x-direction. If set to 'None' it will be chosen automatically.
x_min = None

#Plot upper bound in x-direction. If set to 'None' it will be chosen automatically.
x_max = None

#Plot lower bound in y-direction. If set to 'None' it will be chosen automatically.
y_min = None

#Plot upper bound in y-direction. If set to 'None' it will be chosen automatically.
y_max = None

#Printing interval (i.e. if you don't want to print every step)
print_interval = 1

### END USER INPUT ###





#Importing the basic time data
TDAT = np.loadtxt("Output/time.dat")
nt = TDAT[0]
dt = TDAT[2]


#Checking that the slice direction is allowed and getting nearest grid point
varx = False
vary = False

if variable_x in ['x','X']:
	varx = 'x'
elif variable_x in ['y','Y']:
	varx = 'y'
elif variable_x in ['z','Z']:
	varx = 'z'
elif variable_x in ['vx','VX','Vx','vX']:
	varx = 'vx'
elif variable_x in ['vy','VY','Vy','vY']:
	varx = 'vy'
elif variable_x in ['vz','VZ','Vz','vZ']:
	varx = 'vz'
else:
	print("x-variable is not recognised. Choose \'x\' or \'y\' or \'z\' or \'vx\' or \'vy\' or \'vz\'\nEXITING\n")
	sys.exit()

if variable_y in ['x','X']:
	vary = 'x'
elif variable_y in ['y','Y']:
	vary = 'y'
elif variable_y in ['z','Z']:
	vary = 'z'
elif variable_y in ['vx','VX','Vx','vX']:
	vary = 'vx'
elif variable_y in ['vy','VY','Vy','vY']:
	vary = 'vy'
elif variable_y in ['vz','VZ','Vz','vZ']:
	vary = 'vz'
else:
	print("y-variable is not recognised. Choose \'x\' or \'y\' or \'z\' or \'vx\' or \'vy\' or \'vz\'\nEXITING\n")
	sys.exit()

if(varx == 'z' or vary == 'z'):
	print("Cannot set variable to \'z\' for 2D phase space data\nEXITING\n")
	sys.exit()



#Collecting all of the data files
file_list = glob.glob("Output/" + data_dir + "/*.bin")
nfiles = len(file_list)

#Determing which files correspond to the phase space diagnostic
dlist = []
for i in range(nfiles):
	dstr = re.findall('d\d+',file_list[i])[0]
	dlist.append(int(re.findall('\d+',dstr)[0]))

#Extracting the file list which correspong to the phase space diagnostic
dlist = np.array(dlist)
file_list = np.array(file_list)

file_list = file_list[dlist == diag_num]
nfiles = len(file_list)

if (nfiles == 0):
	print("Phase space diagnostic %d does not exist\nEXITING\n" % (diag_num))
	sys.exit()

#Getting the unique list of time steps from the file names
tlist = []

for i in range(nfiles):
	tstr = re.findall('%d_\d+' % (diag_num),file_list[i])[0]
	tlist.append(int(re.findall('\d+',tstr)[-1]))

tlist = np.sort(tlist)
#tlist = np.array(tlist)


#Getting the diagnostic time increment data
if (len(tlist) > 1):
	ntp = tlist[1] - tlist[0]
else:
	ntp = 1

dtp = ntp*dt


#Getting the start print step. Can be set instead of the start time
ntstart = int(tstart/dt)

#Getting the end print step. Can be set instead of the end time
ntend = int(tend/dt)


#Limiting the start and end time to make sure they fall within the simulation range
if tstart < 0.0:
	tstart = 0.0
	ntstart = 0
	print("Start time cannot be negative, setting to start of simulation:\ntstart = %.4e (s)\n" % (tstart))

if tend < 0.0:
	print("End time cannot be negative\nExiting\n")
	sys.exit()

if tstart > dt*nt:
	print("Start time exceeds simulation end time\nExiting\n")
	sys.exit()

if tend > dt*nt:
	tend = dt*nt
	ntend = tlist[-1]
	print("End time exceeds simulation time, setting to end of simulation:\ntend = %.4e (s)\n" % (tend))



# #Getting the global external grid
XG = np.loadtxt("Output/xgrid.dat")
YG = np.loadtxt("Output/ygrid.dat")
#ZG = np.loadtxt("Output/zgrid.dat")


#Total number of steps to print
NT = tlist[(tlist >= ntstart) & (tlist <= ntend)]
NT = NT[::print_interval]
ntotal = len(NT)
#print(NT)


#Plot settings
label_size = 30
tick_size = 18

#Shifting axis
xshift = 0.06
yshift = 0.02

#Turning off numpy warnings if an empy data set is loaded
warnings.simplefilter("ignore")


if (os.path.isdir(print_dir)):
	os.system('rm -r ' + print_dir)

call(['mkdir',print_dir])


#Working through each time step
for i in range(ntotal):

	fnum = NT[i]

	#Time step progress counter
	print("Step %d of %d" % (i+1,ntotal), end='\r')
	if i+1 == ntotal:
		print("Step %d of %d\n" % (i+1,ntotal))


	DAT = np.fromfile("Output/PhaseDiagnostics/d%d_%08d.bin" % (diag_num,fnum),dtype=np.double)
	npart = int(len(DAT)/5)

	#print(npart)

	#Processing data
	#PSPACE = np.zeros((npart,6))

	x = DAT[0*npart:1*npart]
	y = DAT[1*npart:2*npart]
	#z = DAT[2*npart:3*npart]
	vx = DAT[2*npart:3*npart]
	vy = DAT[3*npart:4*npart]
	vz = DAT[4*npart:5*npart]

	#Setting up the plot based on the x-variable
	if varx == 'x':
		X = x
		xlab = r"$x$ $(m)$"
		xmin = XG[0]
		xmax = XG[-1]
	elif varx == 'y':
		X = y
		xlab = r"$y$ $(m)$"
		xmin = YG[0]
		xmax = YG[-1]
	elif varx == 'z':
		X = z
		xlab = r"$z$ $(m)$"
		xmin = ZG[0]
		xmax = ZG[-1]
	elif varx == 'vx':
		X = vx
		xlab = r"$v_x$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)
	elif varx == 'vy':
		X = vy
		xlab = r"$v_y$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)
	elif varx == 'vz':
		X = vz
		xlab = r"$v_z$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)

	#Setting up the plot based on the y-variable
	if vary == 'x':
		Y = x
		xlab = r"$x$ $(m)$"
		ymin = XG[0]
		ymax = XG[-1]
	elif vary == 'y':
		Y = y
		ylab = r"$y$ $(m)$"
		ymin = YG[0]
		ymax = YG[-1]
	elif vary == 'z':
		Y = z
		ylab = r"$z$ $(m)$"
		ymin = ZG[0]
		ymax = ZG[-1]
	elif vary == 'vx':
		Y = vx
		ylab = r"$v_x$ $(m/s)$"
		ymin = np.min(Y)
		ymax = np.max(Y)
	elif vary == 'vy':
		Y = vy
		ylab = r"$v_y$ $(m/s)$"
		ymin = np.min(Y)
		ymax = np.max(Y)
	elif vary == 'vz':
		Y = vz
		ylab = r"$v_z$ $(m/s)$"
		ymin = np.min(Y)
		ymax = np.max(Y)


	#Setting to user defined min/max if defined
	if (x_min != None):
		xmin = x_min
	if (x_max != None):
		xmax = x_max
	if (y_min != None):
		ymin = y_min
	if (y_max != None):
		ymax = y_max


	#Adjusting the limits for velocity to give some buffer space
	if varx in ['vx','vy','vz']:
		xh = xmax - xmin
		xmax += 0.2*xh
		xmin -= 0.2*xh

	if vary in ['vx','vy','vz']:
		yh = ymax - ymin
		ymax += 0.2*yh
		ymin -= 0.2*yh


	#Plotting
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	plt.scatter(X, Y, c='b', s=20)

	plt.xlabel(xlab, fontsize=label_size)
	plt.ylabel(ylab, fontsize=label_size)

	plt.xlim([xmin,xmax])
	plt.ylim([ymin,ymax])

	plt.tick_params(labelsize=tick_size)


	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	plt.savefig(print_dir + "/%d.png" % i)

	plt.close(fig)
