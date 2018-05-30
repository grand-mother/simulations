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
This script plots the 2D map of the Efield

Usage:  python Efield_2Dmap.py [folder containing the timefresnel-root.dat file] [opt: low frequency cut] [opt: high frequency cut]
"""
    sys.exit(1)
###########################

wkdir = sys.argv[1] # path where the simulation file is

try:
  proj = sys.argv[2]
except:
  proj = 'standart'

try:
  fname = wkdir+'/split/antpos.dat' #ZHAires_output
  print(fname)
  a = np.loadtxt(fname, dtype='float', comments='#')
except:
  fname = wkdir+'/antpos.dat' #ZHAires_output
  a = np.loadtxt(fname, dtype='float', comments='#')

try:
  lowcut = int(sys.argv[3])
  highcut = int(sys.argv[4])
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
for iant in range(0,Nant):
  if lowcut==0 and highcut==0:
    try:
      filename = wkdir+"/split/a"+str(iant)+".trace"
      if not(os.path.exists(filename)):
        filename = wkdir+"/a"+str(iant)+".trace"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti = b[:,0]
      Ampx[iant] = max(b[:,1])-min(b[:,1])
      Ampy[iant] = max(b[:,2])-min(b[:,2])
      Ampz[iant] = max(b[:,3])-min(b[:,3])
    except:
      pass
  else:
    try:
      filename = wkdir+"/split/a"+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      if not(os.path.exists(filename)):
        filename = wkdir+"/a"+str(iant)+"_"+str(lowcut)+"-"+str(highcut)+"MHz.txt"
      b = np.loadtxt(filename, dtype='float', comments='#')
      ti = b[:,0]
      Ampx[iant] = max(b[:,1])-min(b[:,1])
      Ampy[iant] = max(b[:,2])-min(b[:,2])
      Ampz[iant] = max(b[:,3])-min(b[:,3])
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
  pl.scatter(x0,y0,c=Ampx,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ex [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Ex [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Ex_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figx)

  try:
    figx2 = pl.figure(3)
    norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
    #norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
    pl.scatter(x0,y0,c=Ampx,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ex [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ex [$\mu$V/m]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Ex_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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
  pl.scatter(x0,y0,c=Ampy,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ey [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Ey [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Ey_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figy)

  try:
    figy2 = pl.figure(5)
    norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
    #norm = colors.Normalize(vmin=Ampy.min(),vmax=Ampy.max())
    pl.scatter(x0,y0,c=Ampy,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ey [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ey [$\mu$V/m]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Ey_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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
  pl.scatter(x0,y0,c=Ampz,cmap='jet',s=100,norm=norm)
  #pl.title('Amplitude Ez [$\mu$V/m]')
  pl.xlabel(xlbl)
  pl.ylabel(ylbl)
  cbar = pl.colorbar()
  cbar.set_label('Amplitude Ez [$\mu$V/m]')
  #pl.gca().set_aspect(aspect=4) #'equal'
  #pl.axis('equal')
  if lowcut==0 and highcut==0:
    figname = figdir+'/'+showerID+'_ampl_Ez_lin.png'
  else:
    figname = figdir+'/'+showerID+'_ampl_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
  pl.savefig(figname,dpi=500)
  #pl.show()
  #raw_input()
  #pl.close(figz)

  try:
    figz2 = pl.figure(7)
    norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
    #norm = colors.Normalize(vmin=Ampz.min(),vmax=Ampz.max())
    pl.scatter(x0,y0,c=Ampz,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ez [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ez [$\mu$V/m]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    #pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_ampl_Ez_log.png'
    else:
      figname = figdir+'/'+showerID+'_ampl_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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




#################################################################################################################################
# Projection into the shower plane
if proj=='shower':
  sys.path.append('/Users/nrenault/Desktop/GRAND/radiomorphing-master/lib/python/')
  import radiomorphing.frame as frame
  phigeo =0*np.pi/180.  # 182.66#; (ie pointing 2.66 degrees East from full North) # phigeo= 0 from simulations inputfile % In both EVA & Zhaires, North = magnetic North
  thetageo =(180.-27.05)*np.pi/180.  # 152.95*np.pi/180. #27.05*np.pi/180. #; (pointing down)-62.95

  #######################################################################################################################################
  def GetUVW_valentin(pos, core, v, phigeo, bfieldangle):
      relpos = np.array(pos-core)
      inc=bfieldangle

      B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),-np.cos(inc)]) #from oliviers script including phigeo
      B = B/np.linalg.norm(B)
      vxB = np.cross(v,B) #np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
      vxB = vxB/np.linalg.norm(vxB)
      vxvxB = np.cross(v,vxB) #np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
      vxvxB = vxvxB/np.linalg.norm(vxvxB)

      return np.array([np.dot(vxvxB,relpos), np.dot(vxB,relpos), np.dot(v,relpos)]) # vector dot
  #######################################################################################################################################

  azim = 45. #45. #45. #90. #90. #135. #135.
  zen = 100. #100. #105. #100. #100. #100. #105.
  azimr = np.deg2rad(azim)
  zenr = np.deg2rad(zen)
  pos = np.array([x0,y0,z0])
  '''
  R = frame.get_rotation(zenr,azimr,phigeo,thetageo)
  Efield = np.transpose(np.reshape(np.array([Ampx,Ampy,Ampz]),(3,np.shape(np.array([Ampx,Ampy,Ampz]))[1])))
  EshowerA = np.dot(Efield,R)
  Ampx = EshowerA[:,0]#**2
  Ampy = EshowerA[:,1]#**2
  Ampz = EshowerA[:,2]#**2

  GetUVW = frame.UVWGetter(0.,0.,np.mean(z0),zenr,azimr,phigeo,thetageo)

  pos2 = np.zeros((3,len(pos[0,:])))
  for i in np.arange(0,len(pos[0,:])):
    pos2[:,i] = GetUVW(pos[:,i])
  '''
  core = np.reshape(np.array([0.,0.,np.mean(z0)]),(3,1))
  v = -np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), np.cos(zenr)])
  pos2 = GetUVW_valentin(pos, core, v, phigeo, thetageo)

  x0 = pos2[0,:]
  y0 = pos2[1,:]
  z0 = pos2[2,:]
  Ampx=Ampx**2
  Ampy=Ampy**2
  Ampz=Ampz**2


  if DISPLAYX:
    figx1 = pl.figure(2)
    #norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
    norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
    pl.scatter(x0,y0,c=Ampx,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ex [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ex [$\mu$V/m]')
    pl.gca().set_aspect(aspect=4) #'equal'
    pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ex_lin.png'
    else:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
    pl.savefig(figname,dpi=500)
    #raw_input()
    #pl.close(figx)

    try:
      figx2 = pl.figure(3)
      norm = colors.LogNorm(vmin=1, vmax=Ampx.max())
      #norm = colors.Normalize(vmin=Ampx.min(),vmax=Ampx.max())
      pl.scatter(x0,y0,c=Ampx,cmap='jet',s=100,norm=norm)
      #pl.title('Amplitude Ex [$\mu$V/m]')
      pl.xlabel(xlbl)
      pl.ylabel(ylbl)
      cbar = pl.colorbar()
      cbar.set_label('Amplitude Ex [$\mu$V/m]')
      #pl.gca().set_aspect(aspect=4) #'equal'
      pl.axis('equal')
      if lowcut==0 and highcut==0:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ex_log.png'
      else:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ex_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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
    pl.scatter(x0,y0,c=Ampy,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ey [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ey [$\mu$V/m]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ey_lin.png'
    else:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
    pl.savefig(figname,dpi=500)
    #pl.show()
    #raw_input()
    #pl.close(figy)

    try:
      figy2 = pl.figure(5)
      norm = colors.LogNorm(vmin=1, vmax=Ampy.max())
      #norm = colors.Normalize(vmin=Ampy.min(),vmax=Ampy.max())
      pl.scatter(x0,y0,c=Ampy,cmap='jet',s=100,norm=norm)
      #pl.title('Amplitude Ey [$\mu$V/m]')
      pl.xlabel(xlbl)
      pl.ylabel(ylbl)
      cbar = pl.colorbar()
      cbar.set_label('Amplitude Ey [$\mu$V/m]')
      #pl.gca().set_aspect(aspect=4) #'equal'
      pl.axis('equal')
      if lowcut==0 and highcut==0:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ey_log.png'
      else:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ey_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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
    pl.scatter(x0,y0,c=Ampz,cmap='jet',s=100,norm=norm)
    #pl.title('Amplitude Ez [$\mu$V/m]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    cbar = pl.colorbar()
    cbar.set_label('Amplitude Ez [$\mu$V/m]')
    #pl.gca().set_aspect(aspect=4) #'equal'
    pl.axis('equal')
    if lowcut==0 and highcut==0:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ez_lin.png'
    else:
      figname = figdir+'/'+showerID+'_shower_plane_ampl_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_lin.png"
    pl.savefig(figname,dpi=500)
    #pl.show()
    #raw_input()
    #pl.close(figz)

    try:
      figz2 = pl.figure(7)
      norm = colors.LogNorm(vmin=1, vmax=Ampz.max())
      #norm = colors.Normalize(vmin=Ampz.min(),vmax=Ampz.max())
      pl.scatter(x0,y0,c=Ampz,cmap='jet',s=100,norm=norm)
      #pl.title('Amplitude Ez [$\mu$V/m]')
      pl.xlabel(xlbl)
      pl.ylabel(ylbl)
      cbar = pl.colorbar()
      cbar.set_label('Amplitude Ez [$\mu$V/m]')
      #pl.gca().set_aspect(aspect=4) #'equal'
      pl.axis('equal')
      if lowcut==0 and highcut==0:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ez_log.png'
      else:
        figname = figdir+'/'+showerID+'_shower_plane_ampl_Ez_'+str(lowcut)+"-"+str(highcut)+"MHz_log.png"
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
