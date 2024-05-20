#name: iso_b analysis
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

#INPUTS======================
#filepaths

file = 'tide_NC'
for el in np.linspace(100,1000,10):
	year = el

	if file == 'iso_b':
		filename = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/og_data/GLX.0'

	if file == 'tide_b':
		filename = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/tide_b/og_data/GLX.0'

	if file == 'tide_NC':
		filename = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/tide_NC/og_data/GLX.0'

	if 100 <= year <= 600:
		filepth = filename + '00'
		testname = filepth.split('/')[-2]
		timestep = ([str(int(year))])

	if 700 <= year <= 900:
		filepth = filename + '0'
		testname = filepth.split('/')[-2]
		timestep = ([str(int(year))])

	if year == 1000:
		filepth = filename
		testname = filepth.split('/')[-2]
		timestep = (['1000'])

	#Values
	params = {"font.family":"serif","mathtext.fontset":"stix"}
	matplotlib.rcParams.update(params)

	kpc_cgs= 3.08567758e21
	G_cgs  = 6.67e-8
	Mo_cgs = 1.99e33
	umass_GizToGas = 1 #1e9Mo
	umass = 1.0 #* umass_GizToGas
	udist = 1.0  #kpc
	uvel  = np.sqrt( G_cgs * umass * Mo_cgs / (udist * kpc_cgs) )/1e5
	udens = umass * Mo_cgs / (udist * kpc_cgs)**3.
	utime = np.sqrt(1./(udens * G_cgs))
	sec2myr = 60.*60.*24.*365.*1e6

	#access inputs
	i=0
	filenom = (filepth + timestep[i])
	dno = timestep[i]
	s = pynbody.load(filenom)

	s.physical_units()


	t_now =  s.properties['time'].in_units('Myr')
	timestr = str( np.round(float(t_now),1) )
	pynbody.analysis.angmom.faceon(s)

	#BINNING=====================
	r_enc = '15 kpc'
	pp  = pynbody.analysis.profile.Profile(s,max=r_enc,min='0.01 kpc',type='log',nbins=300)
	ppd = pynbody.analysis.profile.Profile(s.dm,max=r_enc,min='0.01 kpc',type='log',nbins=300)
	ppg = pynbody.analysis.profile.Profile(s.gas,max=r_enc,min='0.01 kpc',type='log',nbins=300)
	pps = pynbody.analysis.profile.Profile(s.stars,max=r_enc,min='0.01 kpc',type='log',nbins=300)

	print ('------SETUP-DONE------')

	#PLOTS=======================
	#Galaxy simulation img
	pynbody.plot.image(s.g, width = 100, cmap = 'Blues')
	plt.savefig("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/images/galaxy"+dno+".png", bbox_inches = "tight")
	plt.close()

	# Rotation curve 
	plt.clf()
	plt.plot(pp['rbins'],pp['v_circ'],label='all')
	plt.plot(ppd['rbins'],ppd['v_circ'],label='dark')
	plt.plot(pps['rbins'],pps['v_circ'],label='star')
	plt.plot(ppg['rbins'],ppg['v_circ'],label='gas')
	plt.legend()
	plt.xlim(0,10)
	plt.xlabel('$R$ [kpc]')
	plt.ylabel('$V_{circ}$ [km/s]')
	plt.savefig("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/images/rotation_curve"+dno+".png", bbox_inches = "tight")
	plt.close()

	# Face-on stellar density plot
	vmin=1e0
	vmax=1e4
	plt.clf()
	figS = plt.figure(2)
	axS=figS.add_subplot(1,1,1)
	A1 = sph.image(s.star,qty='rho',width='20 kpc',cmap='bone',units='Msol pc^-2',subplot=axS,vmin=vmin,vmax=vmax,show_cbar = False,resolution=200)
	plt.xlabel('$x \;{\\rm [kpc]}$',fontsize=15)
	plt.ylabel('$y \;{\\rm [kpc]}$',fontsize=15)
	plt.savefig("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/images/stellar_density"+dno+".png", bbox_inches = "tight")
	plt.close()

	# Face-on gas density plot
	vmin=3e-4
	vmax=8e-2
	plt.clf()
	fig=plt.figure(1)
	ax=fig.add_subplot(1,1,1)
	im=sph.image(s.gas,qty='rho',width='30 kpc',cmap='pink',units='g cm^-2',vmin=vmin,vmax=vmax,show_cbar=True,subplot=ax) 
	'''
	#qty = change to any of the derivable,loadable keys. 
	look at the printed statements.
	voids are surrounded by high temperature gas regions.
	stelar ages 'tform'
	only young stars: look at code
	'''
	plt.xlabel('$x \;{\\rm [kpc]}$',fontsize=15)
	plt.ylabel('$y \;{\\rm [kpc]}$',fontsize=15)
	ax.annotate(timestr+'Myr',xy=(0.7,0.9),xycoords='axes fraction',color='white',fontsize=13)
	plt.savefig("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/images/gas_density"+dno+".png", bbox_inches = "tight")
	plt.close()


	#TEXT FILE===================
	Young_Star_names = np.reshape(["snapshot","x", "y", "z", "mass", "tform"],(1,6))
	Young_Star_units = np.reshape(["Myr","kpc", "kpc", "kpc", "$M_{\odot}$", "Myr"],(1,6))
	Young_G_names = np.reshape(["snapshot","x", "y", "z", "mass"],(1,5))
	Young_G_units = np.reshape(["Myr","kpc", "kpc", "kpc", "$M_{\odot}$"],(1,5))
	Table_timeS=np.reshape(np.ones((1,len(s.star['x'])))*int(dno),(len(s.star['x'])))
	Table_timeG=np.reshape(np.ones((1,len(s.gas['x'])))*int(dno),(len(s.gas['x'])))
	print (Table_timeS[0])
	s_forSFR = np.column_stack((Table_timeS,s.star['x'].in_units('kpc'),s.star['y'].in_units('kpc'),s.star['z'].in_units('kpc'),s.star['mass'],s.star['tform'].in_units('Myr')))
	gas_forSFR = np.column_stack((Table_timeG,s.gas['x'].in_units('kpc'),s.gas['y'].in_units('kpc'),s.gas['z'].in_units('kpc'),s.gas['mass']))
	Sout = np.concatenate((Young_Star_names,Young_Star_units,s_forSFR))
	Gas_out = np.concatenate((Young_G_names,Young_G_units,gas_forSFR))
	#print(Sout[1,:])
	np.savetxt("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/og_data/AllStars_"+dno+".txt",Sout,fmt="%s")
	np.savetxt("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/og_data/AllGas_"+dno+".txt",Gas_out,fmt="%s")

	# Stellar velocities for (x,y,z) and (r,phi,z) co-ordinates
	sVel_names = np.reshape(["snapshot","svx", "svy", "svz", "svr", "svphi"],(1,6))
	Vel_units = np.reshape(["Myr","km/s", "km/s", "km/s", "km/s", "km/s"],(1,6))
	Table_timeS=np.reshape(np.ones((1,len(s.star['x'])))*int(dno),(len(s.star['x'])))
	#print Table_timeS[0]
	s_vel = np.column_stack((Table_timeS,s.star['vx'],s.star['vy'],s.star['vz'],s.star['vr'],s.star['vphi']))
	Svout = np.concatenate((sVel_names,Vel_units,s_vel))
	#print(Svout[1,:])
	np.savetxt("/Users/livisilcock/Documents/PROJECTS/VOIDS/data/"+file+"/og_data/AllStarsV_"+dno+".txt",Svout,fmt="%s")





