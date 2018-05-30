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

#pl.ion()
pl.ioff()

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<1 or len(sys.argv)>4):
    print """\
This script plots the 2D map of the difference (absolute and relative) between RadioMorphing computations and ZHAireS simulations either for Efield or voltage.

Usage:  python Comparison_RM2zhaires.py [folder containing the timefresnel-root.dat file] [folder containing the RadioMorphing traces] [Efield or voltage] [opt: low frequency cut] [opt: high frequency cut]
"""
    sys.exit(1)
###########################

VorE = str(sys.argv[3])
if VorE=='Efield':
  prefix = "a"
elif VorE=='voltage':
  prefix = "out_"

try:
  lowcut = int(sys.argv[4])
  highcut = int(sys.argv[5])
except:
  lowcut=0
  highcut=0

############################################################################################################################
## ZHAires part
wkdir = sys.argv[1] # path where the simulation file is
try:
  fname = wkdir+'/split/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')
except:
  fname = wkdir+'/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')

if not(os.path.exists(fname)):
  print "no antpos.dat file"
  print "No antenna within the shower footprint"
  sys.exit(1)

### 
# First load file
x0_zhaires = a[:,0] #X = S->N
y0_zhaires = a[:,1] #Y = E->W
z0_zhaires = a[:,2]  
Nant_zhaires = np.size(x0_zhaires)

Ampx_zhaires = np.zeros((Nant_zhaires,1))
Ampy_zhaires = np.zeros((Nant_zhaires,1))
Ampz_zhaires = np.zeros((Nant_zhaires,1))
for iant in range(0,Nant_zhaires):
  if lowcut==0 and highcut==0:
    try:
      filename = wkdir+"/split/"+prefix+str(iant)+".trace"
      if not(os.path.exists(filename)):
        filename = wkdir+"/"+prefix+str(iant)+".trace"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti_zhaires = b[:,0]
      Ampx_zhaires[iant] = max(abs(b[:,1]))
      Ampy_zhaires[iant] = max(abs(b[:,2]))
      Ampz_zhaires[iant] = max(abs(b[:,3]))
    except:
      pass
  else:
    try:
      filename = wkdir+"/split/"+prefix+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      if not(os.path.exists(filename)):
        filename = wkdir+"/"+prefix+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti_zhaires = b[:,0]
      Ampx_zhaires[iant] = max(abs(b[:,1]))
      Ampy_zhaires[iant] = max(abs(b[:,2]))
      Ampz_zhaires[iant] = max(abs(b[:,3]))
    except:
      pass

# Get amplitude min/max
maxEx_zhaires = np.amax(Ampx_zhaires)
maxEy_zhaires = np.amax(Ampy_zhaires)
maxEz_zhaires = np.amax(Ampz_zhaires)
minEx_zhaires = np.amin(Ampx_zhaires)
minEy_zhaires = np.amin(Ampy_zhaires)
minEz_zhaires = np.amin(Ampz_zhaires)

############################################################################################################################
## RadioMorphing part
wkdir = sys.argv[2] # path where the simulation file is

try:
  fname = wkdir+'/split/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')
except:
  fname = wkdir+'/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')

if not(os.path.exists(fname)):
  print "no antpos.dat file"
  print "No antenna within the shower footprint"
  sys.exit(1)

### 
# First load file
x0_RM = a[:,0] #X = S->N
y0_RM = a[:,1] #Y = E->W
z0_RM = a[:,2]  
Nant_RM = np.size(x0_RM)

Ampx_RM = np.zeros((Nant_RM,1))
Ampy_RM = np.zeros((Nant_RM,1))
Ampz_RM = np.zeros((Nant_RM,1))
for iant in range(0,Nant_RM):
  if lowcut==0 and highcut==0:
    try:
      filename = wkdir+"/split/"+prefix+str(iant)+".trace"
      if not(os.path.exists(filename)):
        filename = wkdir+"/"+prefix+str(iant)+".trace"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti_RM = b[:,0]
      Ampx_RM[iant] = max(abs(b[:,1]))
      Ampy_RM[iant] = max(abs(b[:,2]))
      Ampz_RM[iant] = max(abs(b[:,3]))
    except:
      pass
  else:
    try:
      filename = wkdir+"/split/"+prefix+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      if not(os.path.exists(filename)):
        filename = wkdir+"/"+prefix+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti_RM = b[:,0]
      Ampx_RM[iant] = max(abs(b[:,1]))
      Ampy_RM[iant] = max(abs(b[:,2]))
      Ampz_RM[iant] = max(abs(b[:,3]))
    except:
      pass

# Get amplitude min/max
maxEx_RM = np.amax(Ampx_RM)
maxEy_RM = np.amax(Ampy_RM)
maxEz_RM = np.amax(Ampz_RM)
minEx_RM = np.amin(Ampx_RM)
minEy_RM = np.amin(Ampy_RM)
minEz_RM = np.amin(Ampz_RM)

############################################################################################################################
## Computation part
DAmpx = Ampx_RM-Ampx_zhaires
DAmpy = Ampy_RM-Ampy_zhaires
DAmpz = Ampz_RM-Ampz_zhaires
DAmpx_rel = (Ampx_RM-Ampx_zhaires)/Ampx_zhaires
DAmpy_rel = (Ampy_RM-Ampy_zhaires)/Ampy_zhaires
DAmpz_rel = (Ampz_RM-Ampz_zhaires)/Ampz_zhaires

# Get delta_amplitude min/max
maxDAmpx = np.amax(DAmpx)
maxDAmpy = np.amax(DAmpy)
maxDAmpz = np.amax(DAmpz)
maxDAmpx_rel = np.amax(DAmpx_rel)
maxDAmpy_rel = np.amax(DAmpy_rel)
maxDAmpz_rel = np.amax(DAmpz_rel)
minDAmpx = np.amin(DAmpx)
minDAmpy = np.amin(DAmpy)
minDAmpz = np.amin(DAmpz)
minDAmpx_rel = np.amin(DAmpx_rel)
minDAmpy_rel = np.amin(DAmpy_rel)
minDAmpz_rel = np.amin(DAmpz_rel)

maxD = np.amax([maxDAmpx,maxDAmpy,maxDAmpz])
minD = np.amin([minDAmpx,minDAmpy,minDAmpz])

maxD_rel = np.amax([maxDAmpx_rel,maxDAmpy_rel,maxDAmpz_rel])
minD_rel = np.amin([minDAmpx_rel,minDAmpy_rel,minDAmpz_rel])

############################################################################################################################
## Display part
inclin=0.
print "DEBUG::inclination is hard coded (only impact is that for inclin=90deg, x is not replaced by z in plots)."
if inclin==90.:
  x0_zhaires=y0_zhaires
  y0_zhaires=z0_zhaires
  x0_RM=y0_RM
  y0_RM=z0_RM
  xlbl='Y [m]'
  ylbl='Z [m]'
elif inclin==0. or inclin==5. or inclin==10. or inclin==15. or inclin==20. or inclin==45.:
  xlbl='X [m]'
  ylbl='Y [m]'
else:
  print 'DEBUG: Problem with the reading of the slope'
  quit()

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
  figx1 = pl.figure(1)
  #norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
  norm = colors.Normalize(vmin=minD,vmax=maxD)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpx,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ex [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ex [$\mu$V/m]')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vx [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_amp_Ex_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figx)

  figx2 = pl.figure(2)
  #norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
  norm = colors.Normalize(vmin=minD_rel,vmax=maxD_rel)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpx_rel,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ex [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ex relative')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vx relative')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ex_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figx)

if DISPLAYY:
  figy1 = pl.figure(3)
  #norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
  norm = colors.Normalize(vmin=minD,vmax=maxD)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpy,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ey [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ey [$\mu$V/m]')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vy [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_ampl_Ey_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figy)

  figy2 = pl.figure(4)
  #norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
  norm = colors.Normalize(vmin=minD_rel,vmax=maxD_rel)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpy_rel,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ey [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ey relative')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vy relative')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ey_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figy)

if DISPLAYZ:
  figz1 = pl.figure(5)
  #norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
  norm = colors.Normalize(vmin=minD,vmax=maxD)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpz,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ez [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ez [$\mu$V/m]')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vz [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_ampl_Ez_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figz)

  figz2 = pl.figure(6)
  #norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
  norm = colors.Normalize(vmin=minD_rel,vmax=maxD_rel)
  pl.scatter(x0_zhaires,y0_zhaires,c=DAmpz_rel,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ez [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  if VorE=='Efield':
    cbar.set_label('$\Delta$ Ez relative')
  elif VorE=='voltage':
    cbar.set_label('$\Delta$ Vz relative')
  #pl.gca().set_aspect(aspect=4) #'equal'
  pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ez_lin.png'
  else:
    figname = figdir+'/'+showerID+'_delta_ampl_rel_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figz)

if DISPLAYX:
  #raw_input()
  pl.close(figx1)
  pl.close(figy1)
  pl.close(figz1)
  pl.close(figx2)
  pl.close(figy2)
  pl.close(figz2)
