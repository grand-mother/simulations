#import matplotlib
#matplotlib.use('Agg')

import sys
from sys import argv
import glob
import numpy as np
import pylab as pl
import scipy.interpolate as itp
import matplotlib.colors as colors
import os
import re

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<1 or len(sys.argv)>4):
    print """\
This script plots the 2D map of the voltage

Usage:  python voltage_2Dmap.py [folder containing the timefresnel-root.dat file] [opt: low frequency cut] [opt: high frequency cut]
"""
    sys.exit(1)
###########################

wkdir = sys.argv[1] # path where the simulation file is

try:
  fname = wkdir+'/split/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')
except:
  fname = wkdir+'/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')

try:
  lowcut = int(sys.argv[2])
  highcut = int(sys.argv[3])
except:
  lowcut=0
  highcut=0

if not(os.path.exists(fname)):
  print "no antpos.dat file"
  print "No antenna within the shower footprint"
  sys.exit(1)

###
# First load file
x0 = a[:,0] #X = S->N
y0 = a[:,1] #Y = E->W
z0 = a[:,2]
Nant = np.size(x0)

Ampx = np.zeros((Nant))
Ampy = np.zeros((Nant))
Ampz = np.zeros((Nant))
Amptot = np.zeros((Nant))
for iant in range(0,Nant):
  if lowcut==0 and highcut==0:
    try:
      filename = wkdir+"/split/out_"+str(iant)+".txt"
      if not(os.path.exists(filename)):
        filename = wkdir+"/out_"+str(iant)+".txt"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti = b[:,0]
      Ampx[iant] = max(b[:,1])-min(b[:,1]) #max(abs(b[:,1]))
      Ampy[iant] = max(b[:,2])-min(b[:,2]) #max(abs(b[:,2]))
      Ampz[iant] = max(b[:,3])-min(b[:,3]) #max(abs(b[:,3]))
      Amptot[iant] = max(b[:,1]+b[:,2])-min(b[:,1]+b[:,2]) #max(abs(b[:,3]))
    except:
      pass
  else:
    try:
      filename = wkdir+"/split/out_"+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      if not(os.path.exists(filename)):
        filename = wkdir+"/out_"+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti = b[:,0]
      Ampx[iant] = max(b[:,1])-min(b[:,1]) #max(abs(b[:,1]))
      Ampy[iant] = max(b[:,2])-min(b[:,2]) #max(abs(b[:,2]))
      Ampz[iant] = max(b[:,3])-min(b[:,3]) #max(abs(b[:,3]))
    except:
      pass

inclin=0. #90.
print "DEBUG::inclination is hard coded (only impact is that for inclin=90deg, x is not replaced by z in plots)."
if inclin==90.:
  x0=y0
  y0=z0
  xlbl='Y [m]'
  ylbl='Z [m]'
elif inclin==0. or inclin==5. or inclin==10. or inclin==15. or inclin==20. or inclin==45.:
  xlbl='X [m]'
  ylbl='Y [m]'
else:
  print 'DEBUG: Problem with the reading of the slope'
  quit()

# Get amplitude min/max
maxEx = np.amax(Ampx)
maxEy = np.amax(Ampy)
maxEz = np.amax(Ampz)
minEx = np.amin(Ampx)
minEy = np.amin(Ampy)
minEz = np.amin(Ampz)

# Antenna array 2D-plots
DISPLAYX = 1
DISPLAYY = DISPLAYX
DISPLAYZ = DISPLAYX
showerID = os.path.basename(os.path.dirname(wkdir))
rootdir = os.path.dirname(os.path.dirname(wkdir))
config = os.path.basename(os.path.dirname(wkdir))
figdir1 = rootdir+'/fig/'
if not(os.path.isdir(figdir1)):
  os.mkdir(figdir1)
figdir=figdir1+'/'+config+'/'
if not(os.path.isdir(figdir)):
  os.mkdir(figdir)

if DISPLAYX:
  figx1 = pl.figure(2)
  #norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
  norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
  pl.scatter(x0,y0,c=Ampx,cmap='jet',s=10,norm=norm)
  #pl.title('Amplitude Ex [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Vx [$\mu$V]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Vx_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Vx_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figx)

  try:
    figx2 = pl.figure(3)
    norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
    #norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
    pl.scatter(x0,y0,c=Ampx,cmap='jet',s=10,norm=norm)
    #pl.title('Amplitude Ex [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Vx [$\mu$V]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Vx_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Vx_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
    pl.savefig(figname,dpi=500)
    #pl.show()
    #raw_input()
    #pl.close(figx)
  except:
    pass

if DISPLAYY:
  figy1 = pl.figure(4)
  #norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
  norm = colors.Normalize(vmin=Ampy.min(),vmax=Ampy.max())
  pl.scatter(x0,y0,c=Ampy,cmap='jet',s=10,norm=norm)
  #pl.title('Amplitude Ey [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Vy [$\mu$V]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Vy_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Vy_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figy)

  try:
    figy2 = pl.figure(5)
    norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
    #norm = colors.Normalize(vmin=Ampy.min(),vmax=Ampy.max())
    pl.scatter(x0,y0,c=Ampy,cmap='jet',s=10,norm=norm)
    #pl.title('Amplitude Ey [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Vy [$\mu$V]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Vy_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Vy_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
    pl.savefig(figname,dpi=500)
    #pl.show()
    #raw_input()
    #pl.close(figy)
  except:
    pass

if DISPLAYZ:
  figz1 = pl.figure(6)
  #norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
  norm = colors.Normalize(vmin=Ampz.min(),vmax=Ampz.max())
  pl.scatter(x0,y0,c=Ampz,cmap='jet',s=10,norm=norm)
  #pl.title('Amplitude Ez [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Vz [$\mu$V]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Vz_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Vz_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figz)

  try:
    figz2 = pl.figure(7)
    norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
    #norm = colors.Normalize(vmin=Ampz.min(),vmax=Ampz.max())
    pl.scatter(x0,y0,c=Ampz,cmap='jet',s=10,norm=norm)
    #pl.title('Amplitude Ez [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Vz [$\mu$V]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Vz_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Vz_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
    pl.savefig(figname,dpi=500)
    #pl.show()
    #raw_input()
    #pl.close(figz)
  except:
    pass

if DISPLAYX:
  #raw_input()
  pl.close(figx1)
  pl.close(figy1)
  pl.close(figz1)
  pl.close(figx2)
  pl.close(figy2)
  pl.close(figz2)
