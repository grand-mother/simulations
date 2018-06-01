
# this python script should mimic the functionaliyt of processZhairesShowers.m
# this means it reads in the original ZHAires_output timefresnel_root.dat and splits it into single antenna files a_#.dat. It creates also antpos.dat with all antenna positions
# for filtering and hilbert envelope, please use the next script

#import matplotlib
#matplotlib.use('Agg')

import sys
from sys import argv
import glob
import numpy as np
import pylab as pl
import scipy.interpolate as itp
import os
import re

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<2 or len(sys.argv)>2):
    print """\
This script split the timefresnel-root.dat file

Usage:  python splitZhairesFields_all.py [folder containing the timefresnel-root.dat file]
"""
    sys.exit(1)
###########################

wkdir = sys.argv[1] #path where the simulation file is
fname = wkdir+'timefresnel-root.dat' #ZHAires_output
if not(os.path.exists(fname)):
  print "no timefresnel-root.dat file"
  print "No antenna within the shower footprint"
  sys.exit(1)
if not(os.path.isdir(wkdir+'/split/')):
  os.mkdir(wkdir+'/split/')

print 'file ' +fname+ ' gets split up into single antenna file'

# Get slope from MasterIndex file
master_file=glob.glob(wkdir+'MasterIndex')
if len(master_file)==1:
  try:
      try:
          data_masterfile = file(master_file[0])
          tmp = data_masterfile.readlines()[1]
          alpha=float(tmp.split(' ',-1)[4])
      except:
          data_masterfile = file(master_file[0])
          tmp = data_masterfile.readlines()[2]
          alpha=float(tmp.split(' ',-1)[4])
  except:
      exit()
      alpha = map(int,re.findall('\d+',os.path.basename(os.path.dirname(os.path.dirname(wkdir))).split('_')[2]))[0]
else:
  print 'no MasterIndex file. Setting slope to 0'
  alpha = 0.
#alpha = -alpha

###
# Get ground altitude from input file
steerfile_sim=glob.glob(wkdir+'*.inp')
if len(steerfile_sim)==1:
  datafile = file(steerfile_sim[0])
  for line in datafile:
    if 'GroundAltitude' in line:
      ground_alt = float(line.split(' ',-1)[1])
      break
    else:
      ground_alt = 0.
else:
  ground_alt = 0.

###
# Get injection height from input file
if len(steerfile_sim)==1:
  datafile = file(steerfile_sim[0])
  for line in datafile:
    if 'RASPASSHeight' in line:
      injh = float(line.split(' ',-1)[2])
      break
    else:
      injh = 2000.
else:
  injh = 2000.

# Get ground altitude from input file
if len(steerfile_sim)==1:
  datafile = file(steerfile_sim[0])
  for line in datafile:
    if 'PrimaryAzimAngle' in line:
      azim = float(line.split(' ',-1)[1])+180.
      if azim>360.:
          azim = azim-360.
      break
    else:
      azim = 0.
else:
  azim = 0.
beta = azim - 180.

print 'alpha =',alpha,' beta =',beta

# First load file
a = np.loadtxt(fname, dtype='float', comments='#')
shId = a[:,0]
antId =  a[:,1]
x =   a[:,2]  #X = S->N
y =   a[:,3]  #Y = E->W
z =   100e3-a[:,4] #a[:,4]
t =   a[:,5]  #ns
Ex =  a[:,11]  # V/m
Ey =  a[:,12]
print 'Ez and zpos flipped to correct for coordinate system'
Ez = -1.* a[:,13]

# Z origin is at ground (not sea level) and increasing upward
#for i in range(0, len(z)):
#    z[i]= injh -z[i]
#print 'Flipping z axis...'

# file for antenna positions
file_antpos= wkdir+"/split/antpos.dat"
FILE2 = open(file_antpos, "w" )

Nantennas = int(max(antId))
print "Array of ",Nantennas

x0 = np.zeros((Nantennas))
y0 = np.zeros((Nantennas))
z0 = np.zeros((Nantennas))
Ampx = np.zeros((Nantennas))
Ampy = np.zeros((Nantennas))
Ampz = np.zeros((Nantennas))

# Now split antenna data
print 'DEBUG:::Warning: beta is hard coded.'
for i in range(0, Nantennas):
  sel = np.where(antId == i+1)[0]
  x0[i] = x[sel[0]]
  y0[i] = y[sel[0]]
  z0[i] = z[sel[0]]

for i in range(0, Nantennas):
  sel = np.where(antId == i+1)[0]
  try:
   print >>FILE2,"%.2f	%.2f	%.2f %.2f %.2f" % (x0[i], y0[i],z0[i],alpha,beta)

  except IndexError: # catch the error
    continue

  ti = t[sel] #*1e-3  # [ns]
  Exi = Ex[sel]*1e6   # [muV/m]
  Eyi = Ey[sel]*1e6   # [muV/m]
  Ezi = Ez[sel]*1e6   # [muV/m]
  Ampx[i] = max(abs(Exi))
  Ampy[i] = max(abs(Eyi))
  Ampz[i] = max(abs(Ezi))

  filename = wkdir+"/split/a"+str(i)+".trace"
  alld = np.transpose([ti,Exi,Eyi,Ezi])
  np.savetxt(filename,alld,fmt='%.6e')
  ti0 = ti-ti[0]

  mod_i = i+1 % 10
  if mod_i==0:
    DISPLAY = 0
    if DISPLAY:
      fig = pl.figure(1)
      pl.plot(ti0,Exi,label='Ex (North-South)')
      pl.plot(ti0,Eyi,label='Ey (East-West)')
      pl.plot(ti0,Ezi,label='Ez (Up-Down)')
      pl.xlabel('Time [ns]')
      pl.ylabel('Amplitude [muV/m]')
      pl.grid(True)
      pl.legend(loc='lower right')
      pl.show()
      raw_input()
      pl.close(fig)

FILE2.close()
print 'Timefresnel file splitting finished.'
