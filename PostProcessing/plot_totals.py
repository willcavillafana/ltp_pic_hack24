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
data_dir = "Totals"

#Number of the phase space diagnostic corresponding to the number provided in your LTP-PIC input file
diag_num = 1

#Plot lower bound in t-direction. If set to 'None' it will be chosen automatically.
t_min = None

#Plot upper bound in t-direction. If set to 'None' it will be chosen automatically.
t_max = None




### END USER INPUT ###

#Check if the probe print directory exists. Otherwise create it.
print_dir = "TotalPlots"

if (os.path.isdir(print_dir) == False):
	call(['mkdir',print_dir])



#Plot settings
label_size = 30
tick_size = 18
legend_size = 18

#Shifting axis
xshift = 0.06
yshift = 0.02

#Turning off numpy warnings if an empy data set is loaded
warnings.simplefilter("ignore")




#Plotting total particles
fname = "Output/" + data_dir + "/TotalParticles.dat"
if (os.path.exists(fname)):
	print("Printing Total Particles plot")
	DAT = np.loadtxt(fname)

	T = DAT[:,0]
	DAT = DAT[:,1:]

	if (t_min == None):
		tmin = T[0]
	else:
		tmin = t_min

	if (t_max == None):
		tmax = T[-1]
	else:
		tmax = t_max

	#Plotting
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	plt.plot(T, DAT, linewidth=2)

	plt.xlabel(r"$t$ $(s)$", fontsize=label_size)
	plt.ylabel(r"$Total Particles$", fontsize=label_size)

	plt.xlim([tmin,tmax])
	#plt.ylim([ymin,ymax])
	plt.tick_params(labelsize=tick_size)


	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	plt.savefig(print_dir + "/TotalParticles.png")

	plt.close(fig)

else:
	print("Total Particles data does not exist")







#Plotting total momentum
fname = "Output/" + data_dir + "/TotalMomentum.dat"
if (os.path.exists(fname)):
	print("Printing Total Momentum plot")
	DAT = np.loadtxt(fname)

	#Time data
	T = DAT[:,0]

	#Number of species
	nspecies = (np.shape(DAT)[1] - 1) // 3

	PX = np.zeros(np.shape(DAT)[0])
	PY = np.zeros(np.shape(DAT)[0])
	PZ = np.zeros(np.shape(DAT)[0])

	for i in range(nspecies):
		j = 1 + 3*i
		PX += DAT[:,j]
		PY += DAT[:,j+1]
		PZ += DAT[:,j+2]

	if (t_min == None):
		tmin = T[0]
	else:
		tmin = t_min

	if (t_max == None):
		tmax = T[-1]
	else:
		tmax = t_max


	#Plotting
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	plt.plot(T, PX, linewidth=2)
	plt.plot(T, PY, linewidth=2)
	plt.plot(T, PZ, linewidth=2)

	plt.xlabel(r"$t$ $(s)$", fontsize=label_size)
	plt.ylabel(r"$Total Momentum$ $(kg \cdot m/s)$", fontsize=label_size)

	plt.xlim([tmin,tmax])
	#plt.ylim([ymin,ymax])
	plt.tick_params(labelsize=tick_size)

	plt.legend([r"$P_x$",r"$P_y$",r"$P_z$"],fontsize=legend_size)

	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	plt.savefig(print_dir + "/TotalMomentum.png")

	plt.close(fig)

else:
	print("Total Momentum data does not exist")







#Plotting total momentum
fname = "Output/" + data_dir + "/TotalEnergy.dat"
if (os.path.exists(fname)):
	print("Printing Total Energy plot")
	DAT = np.loadtxt(fname)

	#Time data
	T = DAT[:,0]
	PE = DAT[:,1]
	TE = np.copy(PE)

	#Number of species
	nspecies = np.shape(DAT)[1] - 2

	for i in range(nspecies):
		TE += DAT[:,i+2]

	if (t_min == None):
		tmin = T[0]
	else:
		tmin = t_min

	if (t_max == None):
		tmax = T[-1]
	else:
		tmax = t_max

	labels = ["Total","Potential"]

	#Plotting
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	plt.plot(T, TE, linewidth=2)
	plt.plot(T, PE, linewidth=2)

	for i in range(nspecies):
		plt.plot(T, DAT[:,i+2])
		labels.append("Kinetic %d" % (i))

	plt.xlabel(r"$t$ $(s)$", fontsize=label_size)
	plt.ylabel(r"$Energy$ $(J)$", fontsize=label_size)

	plt.xlim([tmin,tmax])
	#plt.ylim([ymin,ymax])
	plt.tick_params(labelsize=tick_size)

	plt.legend(labels,fontsize=legend_size)

	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+xshift,pos1.y0+yshift,pos1.width,pos1.height]
	ax.set_position(pos2)

	plt.savefig(print_dir + "/TotalEnergy.png")

	plt.close(fig)

else:
	print("Total Energy data does not exist")
