#import matplotlib
#matplotlib.use('Agg')
# Taken from Efield_2Dmap.py

import sys
from sys import argv
import glob
import numpy as np
import pylab as pl
import scipy.interpolate as itp
import matplotlib.colors as colors
import os
#import re
from retro.event import EventIterator, EventLogger
from mpl_toolkits.mplot3d import Axes3D

#pl.ion()
pl.ioff()

###########################
###    Test arguments   ###
###########################
if len(sys.argv)!=2:
    print """\
This script plot the 2D map of the Efield & voltages amplitude pattern

Usage:  python plotAmpPattern.py [folder containing the .json file]
"""
    sys.exit(1)
###########################


def amp2DPlot(zvector,tit):
    pl.figure()
    norm = colors.Normalize(vmin=zvector.min(),vmax=zvector.max())
    pl.scatter(decay_pos[0], decay_pos[1],  c='r', marker='h',s=100)
    pl.plot([decay_pos[0], vend[0]],[decay_pos[1], vend[1]],'r')
    pl.scatter(xants,yants,c=zvector,cmap='jet',s=100,norm=norm)
    pl.axis('equal')
    pl.grid(True)
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    pl.title(tit)
    cbar = pl.colorbar()
    cbar.set_label('peak-peak voltage ($\mu$V)')

wkdir = sys.argv[1] # path where the simulation file is

# First load json infos
json_file =  glob.glob(wkdir+'/*'+'.voltage.json')[0]
for evt in EventIterator(json_file):  # Should only be one
  # Tau decay infos
  tau =  evt["tau_at_decay"]
  decay_pos =  tau[2]
  tau_dir = tau[5]
  zen = tau_dir[0]
  az = tau_dir[1]
  cz = np.cos(zen*np.pi/180)
  sz = np.sin(zen*np.pi/180)
  ca = np.cos(az*np.pi/180)
  sa = np.sin(az*np.pi/180)
  k = np.array([ca*sz,sa*sz,cz])  # Direction vector
  vend = k*50000+decay_pos

  # Ants infos
  ants = np.array(evt["antennas"])
  xants = ants[:,0]
  yants = ants[:,1]
  zants = ants[:,2]
  alpha = ants[:,3]
  beta = ants[:,4]

  # Voltage infos
  [bTrig,antsInd,antsIDs] = checkTrig(event) 
  if bTrig:
  antsin = []
  Ampx=[]
  Ampy=[]
  Ampz=[]
  Ampxy = [];
  v = np.array(evt["voltage"])
  for i in range(np.shape(v)[0]):
    antsin.append(int(v[i,0]))  # Index of antennas with radio simulation
    Ampx.append(float(v[i,1]))  # NS arm
    Ampy.append(float(v[i,2]))  # EW arm
    Ampz.append(float(v[i,3]))  # Vert arm
    Ampxy.append(float(v[i,4]))  # EW arm

  antsin = np.array(antsin[antsInd])
  Ampx = np.array(Ampx[antsInd])
  Ampy = np.array(Ampy[antsInd])
  Ampz = np.array(Ampz[antsInd])
  Ampxy = np.array(Ampxy[antsInd])

print  "Decay at position",decay_pos,"in direction (theta,phi)=",tau_dir
ants = ants[antsin]
xants = xants[antsin]
yants = yants[antsin]
zants = zants[antsin]
alpha = alpha[antsin]
beta = beta[antsin]
imax = np.argmax(Ampx)

print  "Nb of antennas simulated:",np.size(ants)
print  "Nb of antennas trigged (NS,EW,Vert,NS+EW):",np.sum(Ampx>th),np.sum(Ampy>th),np.sum(Ampz>th),np.sum(Ampxy>th/np.sqrt(2))
print "Max:",antsin[imax],ants[imax]

# Antenna array 3D-plots (to check geometry)
fig = pl.figure(1)
ax = fig.add_subplot(111, projection='3d')
slop = alpha
norm = colors.Normalize(vmin=slop.min(),vmax=slop.max())
ax.scatter(xants, yants, zants,c=slop, marker='o',cmap='jet',s=100,norm=norm)
ax.scatter(decay_pos[0], decay_pos[1], decay_pos[2], c='r', marker='h',s=100)
ax.plot([decay_pos[0], vend[0]],[decay_pos[1], vend[1]], [decay_pos[2], vend[2]],'r')
xlbl='SN (m)'
ylbl='EW (m)'
pl.xlabel(xlbl)
pl.ylabel(ylbl)
#ax.zlabel("Altitude [m]")
#ax.axis('equal')
#cbar = ax.colorbar()
#cbar.set_label('Slope alpha [deg]')
fig = pl.figure(2)
ax = fig.add_subplot(111, projection='3d')
slop = alpha
norm = colors.Normalize(vmin=slop.min(),vmax=slop.max())
ax.scatter(xants, yants, alpha,c=slop, marker='o',cmap='jet',s=100,norm=norm)
pl.xlabel(xlbl)
pl.ylabel(ylbl)
fig = pl.figure(3)
ax = fig.add_subplot(111, projection='3d')
slop = beta
norm = colors.Normalize(vmin=slop.min(),vmax=slop.max())
ax.scatter(xants, yants, beta,c=slop, marker='o',cmap='jet',s=100,norm=norm)
pl.xlabel(xlbl)
pl.ylabel(ylbl)

# Amplitude patterns @ ground
amp2DPlot(Ampx,'Vx (NS)')
amp2DPlot(Ampy,'Vy (EW)')
amp2DPlot(Ampz,'Vz (Vert)')
amp2DPlot(Ampxy,'Vx+Vy')
pl.show()
