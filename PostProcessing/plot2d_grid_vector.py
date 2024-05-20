import math
import numpy as np
import sys
import glob
import re
import os
from subprocess import call
import matplotlib.pyplot as plt


#Setting the working directory to the top level code folder
script_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = script_dir.split('/')

if script_dir[-1] == "PostProcessing":
	script_dir = script_dir[:-1]

work_dir = '/'.join(script_dir)
os.chdir(work_dir)



### START USER INPUT ###

#Directory in Output/ to search for data. Pick from ['ElectricField','Current']
data_dir = "Current"

#Directory to create for saving time dependent plots. Time averaged plots are saved in the top directory
print_dir = data_dir + "_plots"

#Start time (s)
tstart = 0.0

#End time (s)
tend = 1.0

#Printing interval (i.e. if you don't want to print every step)
print_interval = 1

#Time average. If True, prints a final time-averaged plot of the data
time_average = True

#Contour or plot lower bound. Set to None to fit bounds of the data at each time step
zmin = None

#Contour or plot upper bound. Set to None to fit bounds of the data at each time step
zmax = None

#Number of contour levels
nlevels = 40

#Streamline density
streamline_density = 1

### END USER INPUT ###







#Importing the time data
TDAT = np.loadtxt("Output/time.dat")
nt = TDAT[0]
#ntp = TDAT[1]
dt = TDAT[2]
#dtp = ntp*dt



#Getting the global external grid
XG = np.loadtxt("Output/xgrid.dat")
YG = np.loadtxt("Output/ygrid.dat")

nx = len(XG)
ny = len(YG)

dx = XG[1] - XG[0]
dy = YG[1] - YG[0]

#Setting up the global 2D grid
XX, YY = np.meshgrid(XG,YG)




#Checking if there is a periodic direction in the simulation
PER = np.loadtxt("Output/periodic.dat")

if PER[0] == 0:
	print("Periodic in X\n")

if PER[1] == 0:
	print("Periodic in Y\n")




#List of all possible data directories for plotting
DATA_dir = ['ElectricField','Current']

#Setting the two directories to collect the vector coordinates from
if data_dir == 'ElectricField':
	dirx = 'ElectricFieldX'
	diry = 'ElectricFieldY'
elif data_dir == 'Current':
	dirx = 'CurrentX'
	diry = 'CurrentY'
else:
	print("This script is not set up to plot this directory. Please pick from one of:")
	print(DATA_dir)
	print("EXITING\n")
	sys.exit()

#Index of the requested data directory
idir = np.where(np.array(DATA_dir) == data_dir)[0][0]



#Collecting all of the data files for the x-component
file_list_x = glob.glob("Output/" + dirx + "/*.bin")
nfilesx = len(file_list_x)

#Getting the unique list of time steps and Regions for the x-components
tlistx = []

for i in range(nfilesx):
	tlistx.append(int(re.findall('\d+',file_list_x[i])[0]))

tlistx = np.unique(tlistx)



#Collecting all of the data files for the y-component
file_list_y = glob.glob("Output/" + diry + "/*.bin")
nfilesy = len(file_list_y)

#Getting the unique list of time steps and Regions for the y-components
tlisty = []

for i in range(nfilesy):
	tlisty.append(int(re.findall('\d+',file_list_y[i])[0]))

tlisty = np.unique(tlisty)


#Checking if we have equal lists
if np.array_equal(tlistx,tlisty) == False:
	print("List of time steps is not identical between the x and y components. Check data.\nEXITING\n")
	sys.exit()

#Creating one list for both directories
tlist = np.copy(tlistx)


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



#Determining whether to set plot contour or y-axis limits
if zmin == None or zmax == None:
	fixed_levels = False
	print("Plot limits will adjust with each time step\n")
else:
	levels = np.linspace(zmin,zmax,nlevels)
	fixed_levels = True
	print("Plot limits are fixed to min = %.4e and max = %.4e\n" % (zmin,zmax))


#Total number of steps to print
NT = tlist[(tlist >= ntstart) & (tlist <= ntend)]
NT = NT[::print_interval]
ntotal = len(NT)


#Plot settings
label_size = 30
tick_size = 18
cbar_size = 14

#Shifting axis
xshift = 0.02
yshift = 0.02

#y-labels
Ylabels = [r'$\mathbf{E}$ $(V/m)$',r'$\mathbf{J}$ $(A/m^2)$']




### Setting up the directory for plots, check if ok to over-write. Automatically overwrites
if (os.path.isdir(print_dir)):
	os.system('rm -r ' + print_dir)

call(['mkdir',print_dir])

#Arrays for storing time averaged data
DATx_tavg = np.zeros([ny,nx])
DATy_tavg = np.zeros([ny,nx])


#Iterating through the time steps
for i in range(ntotal):

	fnum = NT[i]

	#Time step progress counter
	print("Step %d of %d" % (i+1,ntotal), end='\r')
	if i+1 == ntotal:
		print("Step %d of %d\n" % (i+1,ntotal))


	fnamex = "Output/" + dirx + "/%08d.bin" % (tlist[i])
	fnamey = "Output/" + diry + "/%08d.bin" % (tlist[i])

	DATx = np.fromfile(fnamex,dtype=np.double,count=nx*ny)
	DATy = np.fromfile(fnamey,dtype=np.double,count=nx*ny)

	#Reshaping the data to the grid
	DATx = np.reshape(DATx,(ny,nx))
	DATy = np.reshape(DATy,(ny,nx))

	#Handling periodic edges of domain
	if (PER[0] == 0):
		xedge0 = np.array(DATx[:,0])
		DATx[:,-1] = xedge0
		xedge0 = np.array(DATy[:,0])
		DATy[:,-1] = xedge0

	if (PER[1] == 0):
		yedge0 = np.array(DATx[0,:])
		DATx[-1,:] = yedge0
		yedge0 = np.array(DATy[0,:])
		DATy[-1,:] = yedge0


	#Adding to the time averaged data set
	DATx_tavg += DATx
	DATy_tavg += DATy


	#Getting the vector quantity magnitude
	DATmag = (DATx**2 + DATy**2)**0.5

	#Initializing the figure
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	if fixed_levels == False:
		cs = plt.contourf(XX, YY, DATmag, nlevels)
	if fixed_levels == True:
		cs = plt.contourf(XX, YY, DATmag, levels, extend='both')


	plt.streamplot(XX,YY,DATx,DATy,color='orangered')

	plt.xlabel(r"$x$", fontsize=label_size)
	plt.ylabel(r"$y$", fontsize=label_size)

	plt.xlim([XG[0],XG[-1]])
	plt.ylim([YG[0],YG[-1]])

	plt.tick_params(labelsize=tick_size)

	cb = plt.colorbar(cs)
	cb.set_label(label=Ylabels[idir],size=tick_size)
	cb.ax.tick_params(labelsize=cbar_size)

	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	pos1 = cb.ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	cb.ax.set_position(pos2)


	plt.savefig(print_dir + "/%d.png" % i)
	plt.close(fig)


#Time averaging the data
DATx_tavg /= ntotal
DATy_tavg /= ntotal

#Getting the vector quantity magnitude
DATmag_tavg = (DATx_tavg**2 + DATy_tavg**2)**0.5

if (time_average):

	#Initializing the figure
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	if fixed_levels == False:
		cs = plt.contourf(XX, YY, DATmag_tavg, nlevels)
	if fixed_levels == True:
		cs = plt.contourf(XX, YY, DATmag_tavg, levels, extend='both')

	plt.streamplot(XX,YY,DATx_tavg,DATy_tavg,color='orangered')

	plt.xlabel(r"$x$", fontsize=label_size)
	plt.ylabel(r"$y$", fontsize=label_size)

	plt.xlim([XG[0],XG[-1]])
	plt.ylim([YG[0],YG[-1]])

	plt.tick_params(labelsize=tick_size)

	cb = plt.colorbar(cs)
	cb.set_label(label=Ylabels[idir],size=tick_size)
	cb.ax.tick_params(labelsize=cbar_size)

	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	pos1 = cb.ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	cb.ax.set_position(pos2)

	plt.savefig(print_dir + "_tavg.png")
	plt.close(fig)
