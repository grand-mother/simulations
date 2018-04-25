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
if len(sys.argv)<3:
    print """\
This script plot the 2D map of the Efield & voltages amplitude pattern

Usage:  python plotAmpPattern.py [folder containing the .traces file]  [file extension]
"""
    sys.exit(1)
###########################

wkdir = sys.argv[1] # path where the simulation file is

# First load json infos
try:
  json_file =  glob.glob(wkdir+'*'+'.voltage.json')[0]
  for evt in EventIterator(json_file):
    # Ants infos
    ants = np.array(evt["antennas"])
except:
  try: 
    print wkdir
    ants = np.loadtxt(wkdir+"/antpos.dat")    
  except:
    print "caca"
    
    
xants = ants[:,0]
yants = ants[:,1]
zants = ants[:,2]
alpha = ants[:,3]
beta = ants[:,4]
Nant = np.size(xants)

# Now load traces
Ampx = np.zeros((Nant,1))
Ampy = np.zeros((Nant,1))
Ampz = np.zeros((Nant,1))
for iant in range(0,Nant):
  if sys.argv[2]=="trace":
    filename = wkdir+"/a"+str(iant)+".trace"
  else:
     filename = wkdir+"/out_"+str(iant)+".txt"
  if os.path.isfile(filename):
    b = np.loadtxt(filename, dtype='float', comments='#')
    ti = b[:,0]
    #pl.figure(iant)
    #for i in range(1,4):
    #  pl.plot(ti-min(ti),b[:,i])
    #  pl.grid(True)
    #  pl.xlabel('Time (ns)')
    #  pl.ylabel('Efield ($\mu$V/m)')
    #pl.title('Antenna {0}: ({1},{2})'.format(iant,x0[iant],y0[iant]))
    #pl.show()
    Ampx[iant] = max(b[:,1])-min(b[:,1])
    Ampy[iant] = max(b[:,2])-min(b[:,2])
    Ampz[iant] = max(b[:,3])-min(b[:,3])
    print 'Antenna',iant,': max(Ex,Ey,Ez):',Ampx[iant],Ampy[iant],Ampz[iant]

# Antenna array 2D-plots
DISPLAYX = 1

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xants, yants, zants, c='r', marker='o')
pl.xlabel("SN [m]")
pl.ylabel("EW [m]")

if DISPLAYX:
  figx1 = pl.figure(12)
  #norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
  norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
  pl.scatter(xants,yants,c=Ampx,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ex [$\mu$V/m]')
  pl.xlabel("SN [m]")
  pl.ylabel("EW [m]")
  for i in range(len(xants)):
    pl.annotate(str(i),(xants[i]+100,yants[i]+100))
  cbar = pl.colorbar()
  cbar.set_label('peak-peak Ex [$\mu$V/m]')
  figname = 'Ex_lin.png'
  pl.savefig(figname,dpi=500)

  figy1 = pl.figure(13)
  norm = colors.Normalize(vmin=Ampy.min(),vmax=Ampy.max())
  pl.scatter(xants,yants,c=Ampy,cmap='jet',s=100,norm=norm)
  pl.xlabel("SN [m]")
  pl.ylabel("EW [m]")
  cbar = pl.colorbar()
  cbar.set_label('peak-peak Ey [$\mu$V/m]')
  figname = 'Ey_lin.png'
  pl.savefig(figname,dpi=500)

  figz1 = pl.figure(14)
  norm = colors.Normalize(vmin=Ampz.min(),vmax=Ampz.max())
  pl.scatter(xants,yants,c=Ampz,cmap='jet',s=100,norm=norm)
  pl.xlabel("SN [m]")
  pl.ylabel("EW [m]")
  cbar = pl.colorbar()
  cbar.set_label('peak-peak Ez [$\mu$V/m]')
  figname = 'Ez_lin.png'
  pl.savefig(figname,dpi=500)

if DISPLAYX:
  pl.show()
  raw_input()
  pl.close(figx1)
  pl.close(figy1)
  pl.close(figz1)

