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
print_dir = "VDF" + str(diag_num) + "_plots"

#Start time (s)
tstart = 0.0

#End time (s)
tend = 1.0

#Select the velocity coordinate for the VDF. Can be 'vx' or 'vy' or 'vz'
variable = 'vz'

#Plot lower bound. If set to 'None' it will be chosen automatically.
x_min = None

#Plot upper bound. If set to 'None' it will be chosen automatically.
x_max = None

#Number of histogram bins to use
num_bins = None

#Printing interval (i.e. if you don't want to print every step)
print_interval = 2

### END USER INPUT ###





#Importing the basic time data
TDAT = np.loadtxt("Output/time.dat")
nt = TDAT[0]
dt = TDAT[2]


#Checking that the slice direction is allowed and getting nearest grid point
var = False

if variable in ['vx','VX','Vx','vX']:
	var = 'vx'
elif variable in ['vy','VY','Vy','vY']:
	var = 'vy'
elif variable in ['vz','VZ','Vz','vZ']:
	var = 'vz'
else:
	print("Variable is not recognised. Choose \'vx\' or \'vy\' or \'vz\'\nEXITING\n")
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
# XG = np.loadtxt("Output/xgrid.dat")
# YG = np.loadtxt("Output/ygrid.dat")
# ZG = np.loadtxt("Output/zgrid.dat")


#Total number of steps to print
NT = tlist[(tlist >= ntstart) & (tlist <= ntend)]
NT = NT[::print_interval]
ntotal = len(NT)


#Plot settings
label_size = 30
tick_size = 18

#Shifting axis
xshift = 0.04
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

	#x = DAT[0*npart:1*npart]
	#y = DAT[1*npart:2*npart]
	vx = DAT[2*npart:3*npart]
	vy = DAT[3*npart:4*npart]
	vz = DAT[4*npart:5*npart]

	#Setting up the plot based on the x-variable
	if var == 'vx':
		X = vx
		xlab = r"$v_x$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)
	elif var == 'vy':
		X = vy
		xlab = r"$v_y$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)
	elif var == 'vz':
		X = vz
		xlab = r"$v_z$ $(m/s)$"
		xmin = np.min(X)
		xmax = np.max(X)

	#Setting to user defined min/max if defined
	if (x_min != None):
		xmin = x_min
	if (x_max != None):
		xmax = x_max

	#Adjusting the limits for velocity to give some buffer space
	# xh = xmax - xmin
	# xmax += 0.2*xh
	# xmin -= 0.2*xh

	#Producing histogram data
	if (num_bins == None):
		hist, bin_edges = np.histogram(X,range=(xmin,xmax),density=True)
	else:
		hist, bin_edges = np.histogram(X,bins=num_bins,range=(xmin,xmax),density=True)

	#Scaling the histogram data
	bin_length = bin_edges[1] - bin_edges[0]
	hist *= bin_length

	#Getting the bin centers
	bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])


	#Plotting
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	plt.plot(bin_centers, hist, linewidth=2)

	plt.xlabel(xlab, fontsize=label_size)
	plt.ylabel("Probability Density", fontsize=label_size)

	plt.xlim([xmin,xmax])
	#plt.ylim([ymin,ymax])

	plt.tick_params(labelsize=tick_size)


	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	plt.savefig(print_dir + "/%d.png" % i)

	plt.close(fig)
