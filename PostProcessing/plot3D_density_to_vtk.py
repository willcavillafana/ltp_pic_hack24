import math
import numpy as np
import sys
import glob
import re
import os
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from evtk.hl import imageToVTK


#Directories for data and printing
data_dir = "Density"
print_image_dir = "nvtk"


#Plot settings
label_size = 30
tick_size = 18
cbar_size = 14


#Conditions for periodic BCs
PER = np.loadtxt("Output/periodic.dat")

if PER[0] == 0:
	print("Periodic in X")

if PER[1] == 0:
	print("Periodic in Y")
	
if PER[2] == 0:
	print("Periodic in Z")




#Getting the global external grid
XG = np.loadtxt("Output/xgrid.dat")
YG = np.loadtxt("Output/ygrid.dat")
ZG = np.loadtxt("Output/zgrid.dat")

nx = len(XG)
ny = len(YG)
nz = len(ZG)

dx = XG[1] - XG[0]
dy = YG[1] - YG[0]
dz = ZG[1] - ZG[0]

#Setting up the global 2D grid
XX, YY = np.meshgrid(YG,XG)



#Collecting all of the external region grids
file_list_x = glob.glob("Output/Xgrid/*.dat")
nreg_x = len(file_list_x)

file_list_y = glob.glob("Output/Ygrid/*.dat")
nreg_y = len(file_list_y)

file_list_z = glob.glob("Output/Zgrid/*.dat")
nreg_z = len(file_list_z)

if ((nreg_x != nreg_y) or (nreg_x != nreg_z) or (nreg_y != nreg_z)):
	print("Number of x, y or z region grids is not equal, exiting")
	sys.exit()

nreg = nreg_x


#Ordering the file lists
file_num = np.zeros(nreg)

for i in range(nreg):
        file_num[i] = int(re.findall('\d+',file_list_x[i])[0])

file_num = np.sort(file_num)



#Collecting and storing the grids
Xext = []
Yext = []
Zext = []

for i in range(nreg_x):
	fnamex = "Output/Xgrid/xgrid_r%03d.dat" % file_num[i]
	fnamey = "Output/Ygrid/ygrid_r%03d.dat" % file_num[i]
	fnamez = "Output/Zgrid/zgrid_r%03d.dat" % file_num[i]
	Xext.append(np.loadtxt(fnamex))
	Yext.append(np.loadtxt(fnamey))
	Zext.append(np.loadtxt(fnamez))


#print Xext
#print Yext


#Collecting all of the charge data
file_list = glob.glob("Output/" + data_dir + "/*.dat")
nfiles = len(file_list)

slist = []
stlist = []
rlist = []
rtlist = []
tlist = []
ttlist = []


for i in range(nfiles):
        stlist.append(re.findall('s\d+',file_list[i])[0])
        rtlist.append(re.findall('r\d+',file_list[i])[0])
        ttlist.append(re.findall('_\d+',file_list[i])[0])

for i in range(nfiles):
        slist.append(int(re.findall('\d+',stlist[i])[0]))
        rlist.append(int(re.findall('\d+',rtlist[i])[0]))
        tlist.append(int(re.findall('\d+',ttlist[i])[0]))


slist = np.unique(slist)
rlist = np.unique(rlist)
tlist = np.unique(tlist)


#print slist
#print rlist
#print tlist


#Number of time steps
nt = len(tlist)
nt0 = 198
nt1 = nt



#sys.exit()
for sn in slist:
#for sn in [slist[0]]:

	print("Working on species " + str(sn) + " of " + str(len(slist)))

	print_image_dir_species = print_image_dir + str(sn)

	### Setting up the directory for plots, check if ok to over-write
	if (os.path.isdir(print_image_dir_species)):
	        #os.system('rm -r -I ' + print_image_dir)
        	os.system('rm -r ' + print_image_dir_species)

	call(['mkdir',print_image_dir_species])



	### Plotting charge density
	for i in range(nt0,nt1):


		print("Step %d of %d, species %d" % (i+1,nt,sn))

		#Empty grid for storing values
		RHO = np.zeros([nx,ny,nz])

	

		#Collecting data from the individual regions
		for r in range(nreg):

			X = np.array(Xext[r])
			Y = np.array(Yext[r])
			Z = np.array(Zext[r])
			
			nxl = len(X)
			nyl = len(Y)
			nzl = len(Z)

			XGloc = np.where((XG >= X[0]) & (XG <= X[-1]))[0]
			YGloc = np.where((YG >= Y[0]) & (YG <= Y[-1]))[0]
			ZGloc = np.where((ZG >= Z[0]) & (ZG <= Z[-1]))[0]
		
			#rows, cols = np.meshgrid(XGloc, YGloc)

			#xp, yp, zp = np.meshgrid(XGloc,YGloc,ZGloc)
			
			xp, yp = np.meshgrid(XGloc,YGloc)

			#print rows
			#print cols

			#print(ZGloc)

			#Loading the data
			fname = "Output/" + data_dir + "/s%01d_r%03d_%08d.dat" % (sn,r,tlist[i])
			RHOloc = np.loadtxt(fname)/dx**3

			#print np.shape(RHOloc)

			#Storing the data in the global grid on z-level at a time
			for zp in ZGloc:
				kmin = zp*nyl
				kmax = (zp+1)*nyl
				RHO[xp,yp,zp] = RHO[xp,yp,zp] + RHOloc[kmin:kmax,:]

	

		#continue
		#Handling periodic edges of domain
		if (PER[0] == 0):
			xface0 = RHO[0,:,:]
			xface1 = RHO[-1,:,:]

			RHO[0,:,:] += xface1
			RHO[-1,:,:] += xface0


		if (PER[1] == 0):
			yface0 = np.array(RHO[:,0,:])
			yface1 = np.array(RHO[:,-1,:])

			RHO[:,0,:] += yface1
			RHO[:,-1,:] += yface0
			
		if (PER[2] == 0):
			zface0 = RHO[:,:,0]
			zface1 = RHO[:,:,-1]
			
			RHO[:,:,0] += zface1
			RHO[:,:,-1] += zface0



		#print(RHO)
		#continue
		#sys.exit()

		#Exporting to VTK format		
		imageToVTK(print_image_dir_species + "/./out" + str(i), pointData = {"density" : RHO} )
		
		











