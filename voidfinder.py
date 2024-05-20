#name: To find the voids
#author: Olivia Silcock
#date: April 2024


#IMPORTS=========================================
import pynbody
import pynbody.plot.sph as sph
import pynbody.plot as pp
import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import matplotlib
import os
import pynbody.plot.sph as sph
import scipy.signal
import scipy.interpolate
import re
import struct
import pandas as pd

#INPUTS==========================================
#Values
params = {"font.family":"serif","mathtext.fontset":"stix"}
matplotlib.rcParams.update(params)

kpc_cgs= 3.08567758e21
G_cgs  = 6.67e-8
Mo_cgs = 1.99e33
umass_GizToGas = 1.  #1e9Mo
umass = 1.0 #* umass_GizToGas
udist = 1.0  #kpc
uvel  = np.sqrt( G_cgs * umass * Mo_cgs / (udist * kpc_cgs) )/1e5
udens = umass * Mo_cgs / (udist * kpc_cgs)**3.
utime = np.sqrt(1./(udens * G_cgs))
sec2myr = 60.*60.*24.*365.*1e6

#inputs
filepth = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/GLX.0'
testname = filepth.split('/')[-2]
timestep = (['1000'])

#access inputs
i=0
filenom = (filepth + timestep[i])
dno = timestep[i]
s = pynbody.load(filenom)

#params = glob.glob('*.nml')[0] 
print ('-------PREAMBLE-------')
s.physical_units()
print (s.families())
print (s.loadable_keys())
print (s.properties)
print (s.g)
print (s.gas.loadable_keys())
print (s.s)
print (s.star.loadable_keys())
print (s.dm.loadable_keys())

t_now =  s.properties['time'].in_units('Myr')
timestr = str( np.round(float(t_now),1) )
pynbody.analysis.angmom.faceon(s)

#GAS PLOT====================
#Galaxy simulation img
s.physical_units()

vmin=1e-100
vmax=3e-3

plt.clf()
fig=plt.figure(1)
ax=fig.add_subplot(1,1,1)
im=sph.image(s.gas,qty='rho',width='30 kpc',cmap='Set1',units='g cm^-2', vmin = vmin, vmax = vmax, show_cbar=True,subplot=ax) #proj
plt.savefig('/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/images/gas_density_inverted.png', bbox_inches = 'tight')
plt.close()

#VOIDFINDING
#data - this has the predetermined voids' information
tbl = np.loadtxt('/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/VOIDS.csv', delimiter = ',')

#sphere loop to determine each void's scatterplot and minimum density coordinate
for el in range(0,19,1):
	data = tbl[el,:]
	mask = pynbody.filt.Sphere(data[3], cen = (data[0], data[1], data[2]))

	void = s.gas[mask]
	minimum = min(void['rho'])
	print("Minimum is:", minimum)
	min_idx = np.argmin(void['rho'])
	print("The coordinates are:", void['pos'][min_idx, 0], void['pos'][min_idx, 1], void['pos'][min_idx, 2])

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')


	x = void['pos'][:,0]
	y = void['pos'][:,1]
	z = void['pos'][:,2]

	sc = ax.scatter(x, y, z, c=void['rho'], cmap='summer')

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')


	cbar = plt.colorbar(sc)
	cbar.set_label('Density')
	plt.title('Scatter Plot with Density Colormap')
	plt.savefig('/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/VOIDS_scatterplots/'+str(el+1)+'.png')
	plt.close()


#cuboid loop - WIP
tbl = np.loadtxt('/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/VOIDS_cuboid.csv', delimiter = ',')

for el in range(0,19,1):
	data = tbl[el,:]
	mask = pynbody.filt.Cuboid(data[1], data[2], -1.5, data[3], data[4], 1.5)

	void = s.gas[mask]


	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')


	x = void['pos'][:,0]
	y = void['pos'][:,1]
	z = void['pos'][:,2]

	sc = ax.scatter(x, y, z, c=void['rho'], cmap='summer')

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')


	cbar = plt.colorbar(sc)
	cbar.set_label('Density')
	plt.title('Scatter Plot with Density Colormap')
	plt.savefig('/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/VOIDS_scatterplots_cuboid/'+str(el+1)+'.png')
	plt.close()


#EXEMPLAR PLOT===============
#cuboid mask
mask = pynbody.filt.Cuboid(-7.458, 8.250, -5, -5.932, 9.957, 5)
data = s.gas[mask]

#density mask
density_mask = (data['rho'] < 1e5)
data = data[density_mask]

#axes
x = data['pos'][:,0]
y = data['pos'][:,1]
z = data['pos'][:,2]

#plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

sc = ax.scatter(x, y, z, c=data['rho'], cmap='summer')
cbar = plt.colorbar(sc)
cbar.set_label('Density')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Scatter Plot with Density Colormap')
plt.close()









