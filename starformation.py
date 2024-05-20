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

#INPUTS======================
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
filepth = '/Users/livisilcock/Documents/PROJECTS/VOIDS/data/iso_b/og_data/GLX.0'
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

#BINNING=====================
r_enc = '15 kpc'
pp  = pynbody.analysis.profile.Profile(s,max=r_enc,min='0.01 kpc',type='log',nbins=300)
ppd = pynbody.analysis.profile.Profile(s.dm,max=r_enc,min='0.01 kpc',type='log',nbins=300)
ppg = pynbody.analysis.profile.Profile(s.gas,max=r_enc,min='0.01 kpc',type='log',nbins=300)
pps = pynbody.analysis.profile.Profile(s.stars,max=r_enc,min='0.01 kpc',type='log',nbins=300)

print ('------SETUP-DONE------')


#PLOTS=======================
vmin=3e-4
vmax=8e-2
plt.clf()
fig=plt.figure(1)
ax=fig.add_subplot(1,1,1)
im=sph.image(s.gas,qty='rho',width='30 kpc',cmap='pink',units='g cm^-2', vmin = vmin, vmax = vmax, show_cbar=True,subplot=ax) #proj

plt.xlabel('$x \;{\\rm [kpc]}$',fontsize=15)
plt.ylabel('$y \;{\\rm [kpc]}$',fontsize=15)
ax.annotate(timestr+'Myr',xy=(0.7,0.9),xycoords='axes fraction',color='white',fontsize=13)
plt.close()




#to find star formation in the latest period of [dt_age]Myr i.e. 100Myr
AGE_LIM = pynbody.units.Unit('1000 Myr')
R_LIM   = pynbody.units.Unit('10 kpc')
Z_LIM   = pynbody.units.Unit('1 kpc')
T_DELAY = pynbody.units.Unit('0 Myr')  #Ignore this time-frame of SF in past
t_off = 0
t_form = s.star['tform'].in_units('Myr') - t_off # * utime/sec2myr
x=s.star['x']
y=s.star['y']
m=s.star['mass']
mg=s.gas['mass']
xg=s.gas['x']
yg=s.gas['y']
h=['z']
L  = 10. #length of the box eg. 10kpc
dx = 0.1 #step size in x
dy = 0.1 #step size in y
nx = int(L/dx) #number of steps in x
ny = int(L/dy) #number of steps in y
dt_age     = AGE_LIM.in_units('Myr')
max_height = Z_LIM.in_units('kpc')
Mbin = np.zeros((nx,ny))#mass array
Mgbin = np.zeros((nx,ny))#mass array
xbin = np.linspace(-L,L,nx)
ybin = np.linspace(-L,L,ny)
area = np.zeros((nx,ny)) #area array
sfrbin = np.zeros((nx,ny)) #SFR array
logsfrbin = np.zeros((nx,ny)) #logSFR array
t1 = t_now-T_DELAY.in_units('Myr')-dt_age
t2 = t_now-T_DELAY.in_units('Myr')
selage = (t_form>=t1) & (t_form<=t2)
young_m = s.star['mass'][selage]
young_x  = s.star['x'][selage]
young_y  = s.star['y'][selage]
young_h =  np.abs(s.star['z'][selage])
sfr=young_m/(dt_age*1e6)
i=0
j=0
while i<nx:
    j=0
    while j<ny:
        x1 = xbin[i]
        x2 = x1+dx
        y1 = ybin[j]
        y2 = y1+dy
        selgrid =  (young_x>x1) & (young_x<x2) & (young_y>y1) & (young_y<y2) & (young_h<max_height)
        Mbin[i,j] = np.sum(young_m[selgrid])
        area[i,j]  = dx*dy
        sfrbin[i,j] = Mbin[i,j]/area[i,j]/(dt_age*1e6)
        logsfrbin[i,j] = np.log10(sfrbin[i,j])
        j+=1
    i+=1
ii=0
jj=0
while ii<nx:
    jj=0
    while jj<ny:
        x1 = xbin[ii]
        x2 = x1+dx
        y1 = ybin[jj]
        y2 = y1+dy
        selgrid =  (xg>x1) & (xg<x2) & (yg>y1) & (yg<y2)
        Mgbin[ii,jj] = np.sum(mg[selgrid])
        jj+=1
    ii+=1
sfebin=sfrbin/Mgbin
logsfebin=np.log10(sfebin)
plt.imshow(np.fliplr(np.rot90(logsfrbin,k=3)),origin='lower',extent=(-L,L,-L,L),cmap='jet')
plt.clim(-2.5,2.5)
plt.colorbar(orientation="vertical", label="$\log\Sigma_{\\rm SFR}\;{\\rm [M_\\odot yr^{-1} kpc^{-2}]}$")
ax=plt.gca()
ax.set_facecolor((0,0,0.5,1))
plt.xlabel('x [kpc]', fontsize=14)
plt.ylabel('y [kpc]', fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.show()