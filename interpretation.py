
#IMPORTS
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr
import numpy as np
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


#DATA
file = 'tide_NC'
shape = 'Sphere'

if file == 'iso_b':
	rawdata_files = [500, 700, 900, 1000]
if file == 'tide_b':
	rawdata_files = [100, 200, 400, 500, 700, 800, 900, 1000]
if file == 'tide_NC':
	rawdata_files = [600, 700, 800, 900]


for el in rawdata_files:

	if el <= 600:
		zeros = '000'
	if 700 <= el <= 900:
		zeros = '00'
	if el == 1000:
		zeros = '0'


	year = str(int(el))

	#directories
	data_path = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/' + file
	simulation_path = data_path + '/og_data/GLX.' + zeros
	summary = data_path + '/'+ shape + '/Summaries/' + year + '.csv'
 
	df = pd.read_csv(summary) 

	#VALUES======================
	major_axis_deg = df['major_axis_deg']
	minor_axis_deg = df['minor_axis_deg']
	z_axis_deg = 1.5


	#HISTOGRAM n.o. voids========
	#bins
	custom_bins = np.arange(4, 13)  # Bins from 1 to 10
	custom_bins = np.append(custom_bins, np.inf)  # Add infinity to capture values > 10
	#bin_labels = [f'{i}-{i+1}' for i in range(1, 11)] + ['>10']  # Bin labels as strings

	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['dist'], bins=custom_bins, edgecolor='black', alpha=0.7)

	# Add labels and title
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Number of Voids')
	plt.title('Number of Voids with increasing Radii')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/novoids' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#HISTOGRAM n.o. particles====
	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['num_particles'], bins=10, edgecolor='black', alpha=0.7)

	# Add labels and title
	plt.xlabel('Number of particles')
	plt.ylabel('Number of Voids')
	plt.title('Histogram of Voids\' Particles')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/noparticles' + year + '.png', bbox_inches = 'tight')
	plt.close()


	#HISTOGRAM mass==============
	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['sum_mass'], bins=10, edgecolor='black', alpha=0.7)

	# Add labels and title
	plt.xlabel('Void Mass')
	plt.ylabel('Number of Voids')
	plt.title('Number of Voids with Increasing Mass')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/noparticles' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#HISTOGRAM min density=======
	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['min_density'], bins=8, edgecolor='black', alpha=0.7)

	# Add labels and title
	#plt.xticks(ticks=range(1, 12), labels=bin_labels, rotation=45)  # Rotate labels for better visibility
	plt.xlabel('Minimum Density Value')
	plt.ylabel('Number of Voids')
	plt.title('Number of Voids with Increasing Density')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/volume' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#HISTOGRAM surface area======
	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['surface_area'], bins=10, edgecolor='black', alpha=0.7)

	# Add labels and title
	#plt.xticks(ticks=range(1, 12), labels=bin_labels, rotation=45)  # Rotate labels for better visibility
	plt.xlabel('Surface Area (kpc^2)')
	plt.ylabel('Number of Voids')
	plt.title('Histogram of Voids\' Surface Area')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/surfacearea' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#HISTOGRAM volume============
	# Plotting the histogram with custom bins
	plt.figure(figsize=(8, 6))
	plt.hist(df['volume'], bins=9, edgecolor='black', alpha=0.7)

	# Add labels and title
	#plt.xticks(ticks=range(1, 12), labels=bin_labels, rotation=45)  # Rotate labels for better visibility
	plt.xlabel(r'Volume ($kpc^{3}$)')
	plt.ylabel('Number of Voids')
	plt.title('Number of Voids with Increasing Volume')

	# Save
	plt.savefig(data_path + '/' + shape + '/analysis/histograms/volume' + year + '.png', bbox_inches = 'tight')
	plt.close()


	# #2D HISTOGRAM!!!=========
	# # Generate random 2D data for hexbin plot
	# x = df['volume']
	# y = df['surface_area']
	 
	# # Creating a 2D histogram (hexbin plot)
	# plt.hexbin(x, y, gridsize=30, cmap='Blues')
	 
	# # Adding labels and title
	# plt.xlabel('Volume (kpc^3)')
	# plt.ylabel('Surface Area (kpc^2)')
	# plt.title('2D Histogram (Hexbin Plot)')
	 
	# # Adding colorbar
	# plt.colorbar()

	# # Save
	# plt.savefig(data_path + '/' + shape + '/analysis/histograms/2dhist' + year + '.png', bbox_inches = 'tight')
	# plt.close()


	#SCATTERPLOT spatial vs density
	#x values
	correlation, _ = pearsonr(df['x_centre'], df['min_x'])
	plt.scatter(df['x_centre'], df['min_x'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('x_centre (kpc)')
	plt.ylabel('min_x (kpc)')
	plt.title('Comparing the centre of the void vs the min. density point')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/x_scatterplot' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#y values
	correlation, _ = pearsonr(df['y_centre'], df['min_y'])
	plt.scatter(df['y_centre'], df['min_y'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('y_centre (kpc)')
	plt.ylabel('min_y (kpc)')
	plt.title('Comparing the centre of the void vs the min. density point')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/y_scatterplot' + year + '.png', bbox_inches = 'tight')
	plt.close()

	#z values
	correlation, _ = pearsonr(df['z_centre'], df['min_z'])
	plt.scatter(df['z_centre'], df['min_z'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('z_centre (kpc)')
	plt.ylabel('min_z (kpc)')
	plt.title('Comparing the centre of the void vs the min. density point')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/z_scatterplot' + year + '.png', bbox_inches = 'tight')

	#SCATTERPLOT volume vs density
	correlation, _ = pearsonr(df['volume'], df['min_density'])
	plt.scatter(df['volume'], df['min_density'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('Volume')
	plt.ylabel('Minium Density')
	plt.title('Comparing the volume to the minimum density value')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/volume_vs_density' + year + '.png', bbox_inches = 'tight')	
	plt.close()

	#SCATTERPLOT n.o. particles vs sum mass
	correlation, _ = pearsonr(df['num_particles'], df['sum_mass'])
	plt.scatter(df['num_particles'], df['sum_mass'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('Number of particles')
	plt.ylabel('Mass')
	plt.title('Comparing the number of particles in the Void vs the Void Mass')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/noparticles_vs_mass' + year + '.png', bbox_inches = 'tight')
	plt.close()

	 #SCATTERLOT volume vs surface area
	correlation, _ = pearsonr(df['volume'], df['surface_area'])
	plt.scatter(df['volume'], df['surface_area'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('Volume')
	plt.ylabel('Surface Area')
	plt.title('Comparing the Void Volume to its Surface Area')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/volume_vs_SA_scatterplot' + year + '.png', bbox_inches = 'tight')
	plt.close()

	 #SCATTERLOT volume vs surface area
	correlation, _ = pearsonr(df['dist'], df['sum_mass'])
	plt.scatter(df['dist'], df['sum_mass'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Sum of Mass')
	plt.title('Void Mass vs Radius')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/mass_vs_radius_scatterplot' + year + '.png', bbox_inches = 'tight')
	plt.close()

	 #SCATTERLOT volume vs surface area
	correlation, _ = pearsonr(df['dist'], df['volume'])
	plt.scatter(df['dist'], df['volume'])
	plt.text(0, 0, f"{correlation}")
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Volume')
	plt.title('Void Volume vs Radius')
	plt.savefig(data_path + '/' + shape + '/analysis/scatterplots/volume_vs_radius_scatterplot' + year + '.png', bbox_inches = 'tight')
	plt.close()







