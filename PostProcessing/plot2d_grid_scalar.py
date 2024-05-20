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

#Directory in Output/ to search for data. Pick from ['Density','ChargeDensity','ChargeDensityBlocks','Potential','ElectricFieldX','ElectricFieldY','CurrentX','CurrentY','CurrentZ','Temperature','TemperatureX','TemperatureY','TemperatureZ','nVX','nVY','nVZ','nVX2','nVY2','nVZ2']
data_dir = "Density"

#Directory to create for saving time dependent plots. Time averaged plots are saved in the top directory
print_dir = data_dir + "_plots"

#Start time (s)
tstart = 0.0

#End time (s)
tend = 1.0

#Printing interval (i.e. if you don't want to print every step)
print_interval = 1

#Time average. If True, prints a final time-averaged plot of the data
time_average = False

#Space averaging. If True must set average_direction to 'x' or 'y'
space_average = False
average_direction = 'x' #'y'

#Species list. List of species to print (refers to numbering of input file). Must be a list, i.e. for one species use species_list = [1]
species_list = [1,2]

#Contour or plot lower bound. Set to None to fit bounds of the data at each time step
zmin = None

#Contour or plot upper bound. Set to None to fit bounds of the data at each time step
zmax = None

#Number of contour levels
nlevels = 40

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
DATA_dir = ['Density','ChargeDensity','ChargeDensityBlocks','Potential','ElectricFieldX','ElectricFieldY','CurrentX','CurrentY','CurrentZ','Temperature','TemperatureX','TemperatureY','TemperatureZ','nVX','nVY','nVZ','nVX2','nVY2','nVZ2']

#Checking that the requested directory is in the above list
if len(np.where(np.array(DATA_dir) == data_dir)[0]) == 0:
	print("This script is not set up to plot this directory. Please pick from one of:")
	print(DATA_dir)
	print("EXITING\n")
	sys.exit()

#Index of the requested data directory
idir = np.where(np.array(DATA_dir) == data_dir)[0][0]


#Determining whether the input is species dependent
species = False
if data_dir in ['Density','Temperature','TemperatureX','TemperatureY','TemperatureZ','nVX','nVY','nVZ','nVX2','nVY2','nVZ2']:
	species = True




#Collecting all of the data files
file_list = glob.glob("Output/" + data_dir + "/*.bin")
nfiles = len(file_list)

#Getting the unique list of time steps, Regions and species from the file names
tlist = []
slist = []

for i in range(nfiles):
	if (species):
		tstr = re.findall('_\d+',file_list[i])[0]
		tlist.append(int(re.findall('\d+',tstr)[0]))
		sstr = re.findall('s\d+',file_list[i])[0]
		slist.append(int(re.findall('\d+',sstr)[0]))
	else:
		tlist.append(int(re.findall('\d+',file_list[i])[0]))


slist = np.unique(slist)
tlist = np.unique(tlist)



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



#Limiting the species to those which exist within the simulation.
if (species):
	species_list = np.array(species_list)
	species_list = np.unique(species_list)

	sind = np.in1d(species_list,slist)
	slist = species_list[sind]

	for i in range(len(species_list)):
		if sind[i] == False:
			print("Species %d will not be processed as it does not exist in the simulation\n" % (species_list[i]))
else:
	slist = [0]


#Determining whether we want to compute spatial averages
if space_average == False:
	xavg = False
	yavg = False
else:
	if average_direction == 'x' or average_direction == 'X':
		xavg = True
		yavg = False
		print_dir += '_xavg'
		print("Performing averages along the x-direction\n")
	elif average_direction == 'y' or average_direction == 'Y':
		xavg = False
		yavg = True
		print_dir += '_yavg'
		print("Performing averages along the y-direction\n")
	else:
		print("Average direction %s is not recognised. Choose \'x\' or \'y\'\nEXITING\n" % (average_direction))
		sys.exit()


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


#sys.exit()
#Plot settings
label_size = 30
tick_size = 18
cbar_size = 14

#Shifting axis
xshift = 0.02
yshift = 0.02

#y-labels
Ylabels = [r'$n$ $(1/m^3)$',r'$\rho$ $(C/m^3)$',r'$\rho$ $(C/m^3)$',r'$\phi$ $(V)$',r'$E_{x}$ $(V/m)$',r'$E_{y}$ $(V/m)$',r'$J_{x}$ $(A/m^2)$',r'$J_{y}$ $(A/m^2)$',r'$J_{z}$ $(A/m^2)$',r'$T$ $(eV)$',r'$T_{x}$ $(eV)$',r'$T_{y}$ $(eV)$',r'$T_{z}$ $(eV)$',r'$nV_{x}$ $(1/m^2 \cdot s)$',r'$nV_{y}$ $(1/m^2 \cdot s)$',r'$nV_{z}$ $(1/m^2 \cdot s)$',r'$nV^2_{x}$ $(1/m \cdot s^2)$',r'$nV^2_{y}$ $(1/m \cdot s^2)$',r'$nV^2_{z}$ $(1/m \cdot s^2)$']






#sys.exit()
#Iterating through the species
for sn in slist:

	if (species):
		print("Working on species " + str(sn) + " of " + str(len(slist)))

	### Setting up the directory for plots, check if ok to over-write. Automatically overwrites
	if (species):
		print_dir_species = print_dir + str(sn)
	else:
		print_dir_species = print_dir

	if (os.path.isdir(print_dir_species)):
		os.system('rm -r ' + print_dir_species)

	call(['mkdir',print_dir_species])


	#Arrays for storing time averaged data
	if (xavg):
		DAT_tavg = np.zeros(ny)
	elif (yavg):
		DAT_tavg = np.zeros(nx)
	else:
		DAT_tavg = np.zeros([ny,nx])
		#DAT_tavg = np.zeros([ny-1,nx-1])


	#Iterating through the time steps
	for i in range(ntotal):

		fnum = NT[i]

		#Time step progress counter
		print("Step %d of %d" % (i+1,ntotal), end='\r')
		if i+1 == ntotal:
			print("Step %d of %d\n" % (i+1,ntotal))


		#Loading the data
		if (species):
			fname = "Output/" + data_dir + "/s%01d_%08d.bin" % (sn,fnum)
		else:
			fname = "Output/" + data_dir + "/%08d.bin" % (fnum)

		DAT = np.fromfile(fname,dtype=np.double,count=nx*ny)

		#Reshaping the data to the grid
		DAT = np.reshape(DAT,(ny,nx))


		# #Handling periodic edges of domain
		if (PER[0] == 0):
			xedge0 = np.array(DAT[:,0])
			DAT[:,-1] = xedge0

		if (PER[1] == 0):
			yedge0 = np.array(DAT[0,:])
			DAT[-1,:] = yedge0


		#Computing relevant averages if requested
		if xavg == True:
			DAT = np.mean(DAT,axis=1)
		if yavg == True:
			DAT = np.mean(DAT,axis=0)


		#Adding to the time averaged data set
		DAT_tavg += DAT


		#Initializing the figure
		fig = plt.figure(figsize=(10,8))
		ax = fig.add_subplot(111)

		#Contour plots
		if xavg == False and yavg == False:
			if fixed_levels == False:
				cs = plt.contourf(XX, YY, DAT, nlevels)
			if fixed_levels == True:
				cs = plt.contourf(XX, YY, DAT, levels, extend='both')

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


		#Plotting x-averaged data
		if xavg == True:
			plt.plot(YG, DAT, linewidth=2)

			plt.xlabel(r"$y$", fontsize=label_size)
			plt.ylabel(Ylabels[idir], fontsize=label_size)

			plt.xlim([YG[0],YG[-1]])

			if fixed_levels == True:
				plt.ylim([zmin,zmax])

			plt.tick_params(labelsize=tick_size)

			#Shifting axis slightly
			pos1 = ax.get_position()
			pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
			ax.set_position(pos2)


		#Plotting y-averaged data
		if yavg == True:
			plt.plot(XG, DAT, linewidth=2)

			plt.xlabel(r"$x$", fontsize=label_size)
			plt.ylabel(Ylabels[idir], fontsize=label_size)

			plt.xlim([XG[0],XG[-1]])

			if fixed_levels == True:
				plt.ylim([zmin,zmax])

			plt.tick_params(labelsize=tick_size)

			#Shifting axis slightly
			pos1 = ax.get_position()
			pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
			ax.set_position(pos2)


		plt.savefig(print_dir_species + "/%d.png" % i)
		plt.close(fig)



	#Time averaging the data
	DAT_tavg /= ntotal

	if (time_average):

		#Initializing the figure
		fig = plt.figure(figsize=(10,8))
		ax = fig.add_subplot(111)

		#Contour plots
		if xavg == False and yavg == False:

			if fixed_levels == False:
				cs = plt.contourf(XX, YY, DAT_tavg, nlevels)
			if fixed_levels == True:
				cs = plt.contourf(XX, YY, DAT_tavg, levels, extend='both')

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

			plt.savefig(print_dir_species + "_tavg.png")

			plt.close(fig)


		#Plotting x-averaged data
		if xavg == True:

			plt.plot(YG, DAT_tavg, linewidth=2)

			plt.xlabel(r"$y$", fontsize=label_size)
			plt.ylabel(Ylabels[idir], fontsize=label_size)

			plt.xlim([YG[0],YG[-1]])

			if fixed_levels == True:
				plt.ylim([zmin,zmax])

			plt.tick_params(labelsize=tick_size)

			#Shifting axis slightly
			pos1 = ax.get_position()
			pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
			ax.set_position(pos2)

			plt.savefig(print_dir_species + "_tavg.png")

			plt.close(fig)


		#Plotting y-averaged data
		if yavg == True:

			plt.plot(XG, DAT_tavg, linewidth=2)

			plt.xlabel(r"$x$", fontsize=label_size)
			plt.ylabel(Ylabels[idir], fontsize=label_size)

			plt.xlim([XG[0],XG[-1]])

			if fixed_levels == True:
				plt.ylim([zmin,zmax])

			plt.tick_params(labelsize=tick_size)

			#Shifting axis slightly
			pos1 = ax.get_position()
			pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
			ax.set_position(pos2)

			plt.savefig(print_dir_species + "_tavg.png")

			plt.close(fig)
