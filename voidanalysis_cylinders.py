#name: To analyse the voids
#author: Olivia Silcock
#date: April 2024


#IMPORTS=====================
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
import csv
import math
from matplotlib.patches import Ellipse

#FUNCTIONS===================
def append_to_csv(data, filename):
    file_exists = os.path.exists(filename)
    with open(filename, mode='a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['x_centre', 'y_centre', 
        	'z_centre', 'dist', 'min_density', 'min_x', 'min_y', 'min_z', 
        	'major_axis_deg', 'minor_axis_deg', 'num_particles', 'sum_mass', 
        	'volume', 'surface_area'])
        if not file_exists:  # Write headers only if the file is newly created
            writer.writeheader()
        writer.writerow(data)


def reorder_csv(filename):
    # Read the existing data from the CSV file into a list of dictionaries
    data = []
    with open(filename, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            data.append(row)

    # Sort the data based on the 'dist' (distance) key in each dictionary
    sorted_data = sorted(data, key=lambda x: float(x['dist']))

    # Write the sorted data back to the CSV file
    with open(filename, mode='w', newline='') as file:
        fieldnames = ['x_centre', 'y_centre', 'z_centre', 'dist', 'min_density', 
        'min_x', 'min_y', 'min_z', 'major_axis_deg', 'minor_axis_deg',
        'num_particles', 'sum_mass', 'volume', 'surface_area']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted_data:
            writer.writerow(row)

#CONSTANTS===================
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

#filenames
file = 'tide_b'

if file == 'iso_b':
	rawdata_files = [500, 700, 900, 1000]
if file == 'tide_b':
	rawdata_files = [100, 200, 400, 500, 700, 800, 900, 1000]
if file == 'tide_NC':
	rawdata_files = [600, 700, 800, 900]


for el in rawdata_files:
	year = el

	if 100 <= year <= 600:
		zeros = '000'

	if 700 <= year <= 900:
		zeros = '00'


	if year == 1000:
		zeros = '0'

	year = str(int(year))

	data_path = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/' + file
	simulation_path = data_path + '/og_data/GLX.' + zeros
	csv_path = data_path + '/Cylinder/Summaries/' + year + '.csv'
	tbl = np.loadtxt(data_path + '/Rawdata/' + year + '.csv', delimiter = ',' )		

	#DATA========================
	testname = simulation_path.split('/')[-2]
	timestep = ([year])

	i=0
	filenom = (simulation_path + timestep[i])
	dno = timestep[i]
	s = pynbody.load(filenom)

	t_now =  s.properties['time'].in_units('Myr')
	timestr = str( np.round(float(t_now),1) )
	pynbody.analysis.angmom.faceon(s)
	s.physical_units()


	#VOID_ANALYSIS===============
	#sphere loop to determine each void's scatterplot and minimum density coordinate
	for el in range(0,len(tbl),1):
		#access the relevant void data
		data = tbl[el,:]

		a = data[5]/2
		b = data[6]/2
		c = 1.5

		x_centre = data[9]
		y_centre = data[10]
		z_centre = data[11]

		#define spherical mask
		mask = pynbody.filt.EllipticalDisc(a, b, c, cen = (x_centre, y_centre, z_centre))
		void = s.gas[mask]

		#minimum density point
		if len(void) == 0:
			minimum = 0
			min_idx = [x_centre, y_centre, z_centre]
		else:
			minimum = min(void['rho'])
			indexing = np.argmin(void['rho'])
			min_idx = void['pos'][indexing]



		kwargs = {
		'x_centre': x_centre,
		'y_centre': y_centre,
		'z_centre': z_centre,
		'dist': math.sqrt((data[0])**2 + (data[1])**2 + (data[2])**2),
		'min_density': minimum,
		'min_x': min_idx[0],
		'min_y': min_idx[1],
		'min_z': min_idx[2],
		'major_axis_deg': a,
		'minor_axis_deg': b,
		'num_particles': len(void),
		'sum_mass': np.sum(void['mass']),
		'volume': np.pi*a*b*3,
		'surface_area': 2*np.pi*a*b + np.pi* 3 * (a + b)
		}

		append_to_csv(kwargs, csv_path)
		reorder_csv(csv_path)


		#PLOTS===================
		#plot - scatterplot
		fig = plt.figure()
		ax = fig.add_subplot(projection = '3d')
		sc = ax.scatter(void['pos'][:,0], void['pos'][:,1], void['pos'][:,2], c=void['rho'], cmap = 'summer')

		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')

		cbar = plt.colorbar(sc)
		cbar.set_label('Density')
		plt.title('Scatter Plot with Density Colormap')
		plt.savefig(data_path + '/Cylinder/scatterplots/' + year  + '/' + str(el) + '.png',)
		plt.close()

		#plot - density
		x = np.linspace(1, 100, 1000)
		y = []

		for number in x:
			major_axis = a/number
			minor_axis = b/number
			z_axis = 1.5/number

			profiling_mask = pynbody.filt.EllipticalDisc(major_axis, minor_axis,
				z_axis, cen = (x_centre, y_centre, z_centre))
			profiling_void = void[profiling_mask]
			density =len(profiling_void)/((4/3)*np.pi*major_axis*minor_axis*1.5)
			y.append(density)

		plt.scatter(x, y, s = 1000, color = 'r', marker = 'x')
		plt.xlabel('Radius(kpc)')
		plt.ylabel('Density(particles/$cm^{3}$)')
		plt.title('Radial Density Profile')
		plt.savefig(data_path + '/Cylinder/radial_density/'+ year+'/' +str(el)+'.png')
		plt.close()
		y.clear()



