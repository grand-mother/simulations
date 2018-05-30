'''
    Welcome to RADIO MORPHING
    start the run with: python example.py
'''

#!/usr/bin/env python
from os.path import split, join, realpath
import os
import sys
import numpy as np
root_dir = realpath(join(split(__file__)[0], ".."))
sys.path.append(join(root_dir, "lib", "python"))
sys.path.append('/Users/nrenault/Desktop/GRAND/retro-master/lib/python/')
sys.path.append('/Users/nrenault/Desktop/GRAND/radiomorphing-master/lib/python/')
import radiomorphing

wkdir = '/Users/nrenault/Desktop/GRAND/scripts_clean/'
GRAND_dir = "/Users/nrenault/Desktop/GRAND/"
RM_dir = join(GRAND_dir,"RadioMorphing")
sim_dir = join(RM_dir, "ref","GrandEventADetailed2") # folder containing your refernence shower simulations
GdAlt = 1500 #m above sea level

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<2 or len(sys.argv)>2):
    print """\
Create the input parameter lists for RM and run the RadioMorphing simulations.

Usage:  python run_RM.py [directory containing the inp zhaires files] [path to the txt file with the shower list]
"""
    sys.exit(1)
###########################

#===========================================================================================================
def inputfromtxt(input_file_path):
#===========================================================================================================
    particule = ['eta','pi+','pi-','pi0','Proton','p','proton','gamma','Gamma','electron','Electron','e-','K+','K-','K0L','K0S','K*+'
    ,'muon+','muon-','Muon+','Muon-','mu+','mu-','tau+','tau-','nu(t)','Positron','positron','e+']

    datafile = file(input_file_path)

    for line in datafile:
        if 'PrimaryZenAngle' in line:
            zen=float(line.split(' ',-1)[1])
            zen = 180-zen  #conversion to GRAND convention i.e. pointing towards antenna/propagtion direction
        if 'PrimaryAzimAngle' in line:
            azim = float(line.split(' ',-1)[1])+180 #conversion to GRAND convention i.e. pointing towards antenna/propagtion direction
            if azim>=360:
                azim= azim-360
        if 'RASPASSHeight' in line:
            injh = float(line.split(' ',-1)[2])
        if 'PrimaryEnergy' in line:
            energy = float(line.split(' ',-1)[1])
        if 'PrimaryParticle' in line:
            primarytype = str(line.split(' ',-1)[1])
            if primarytype[-1]=='\n':
                primarytype=primarytype[0:-1]
        if 'AddSpecialParticle      RASPASSMulti' in line:
            RASPASSMulti_line = line

    try:
        injh
    except NameError:
        injh = 100000. #Case of a cosmic for which no injection height is defined in the input file and is then set to 100 km by ZHAireS
    try:
        energy
    except NameError:
        print 'No primary energy found in the ZHAireS input text file.'
        exit()
    try:
        primarytype
    except NameError:
        primarytype = None

    #energy = energy *1e-18

    if primarytype=='RASPASSMulti':
        tmp = RASPASSMulti_line.split(' ',-1)
        if tmp[-1][-1]=='\n':
            tmp[-1]=tmp[-1][0:-1]
        prod = [x for x in particule if x in set(tmp)]
        ind_prod = np.array([tmp.index(x) for x in prod],dtype=int)
        Wprod = [float(tmp[ind]) for ind in ind_prod+1]
        primarytype = prod[np.argmax(Wprod)]

    if primarytype=='Proton' or primarytype=='K+' or primarytype=='K-' or primarytype=='K0L' or primarytype=='K0S' or primarytype=='K*+':
        primarytype='proton'
    elif primarytype=='gamma' or primarytype=='Gamma' or primarytype=='Electron':
        primarytype='electron'
    elif primarytype=='pi0' or primarytype=='pi-' or primarytype=='pi+':
        primarytype='pion'

    return zen,azim,energy,injh,primarytype

#===========================================================================================================
def create_antpos(out_dir,datafile,inclin):
#===========================================================================================================
  file_antpos= join(out_dir,"antpos.dat")

  if inclin==None:
    FILE2 = open(file_antpos, "w" )
    for line in datafile:
      if 'AddAntenna' in line:
        print >>FILE2,"%.2f %.2f %.2f" % (float(line.split(' ',-1)[1]),float(line.split(' ',-1)[2]),float(line.split(' ',-1)[3])+GdAlt)
  else:
    os.remove(file_antpos)
    FILE2 = open(file_antpos, "w" )
    for line in datafile:
      if 'AddAntenna' in line:
        print >>FILE2,"%.2f %.2f %.2f %.2f %.2f" % (float(line.split(' ',-1)[1]),float(line.split(' ',-1)[2]),float(line.split(' ',-1)[3])+GdAlt,inclin[0],inclin[1])
  return

########################################################################################################################################
########################################################################################################################################
if __name__ == '__main__':

  # Settings of the radiomorphing
  inp_dir = sys.argv[1] #directory where are located the ZHAireS input files
  shower_list = np.loadtxt(sys.argv[2])
  Nsh = np.size(shower_list)
  slope = np.array([-10.,0.],dtype=float)
  task = 'EE1E+10_D40000_alpha10_height3000_sep500'


  for ish in range(0,Nsh):
    shower_name = str(int(shower_list[ish]))
    print shower_name
    inp_file = inp_dir+"D40000m-Z10deg-"+shower_name+".inp"
    try :
        datafile = file(inp_file)

        out_dir = join(RM_dir, "simus_RM",task,shower_name) # folder which will contain radio morphed traces afterwards
        if not(os.path.isdir(out_dir)):
            os.mkdir(out_dir)

        create_antpos(out_dir,datafile,None)

        antennas = join(out_dir, "antpos.dat") # list of antenna positions you would like to simulate, stored in out_dir in the best case
        zenith_sim,azimuth_sim,energy,injection_height,primary = inputfromtxt(inp_file)
        energy = energy/1e18 #convert energy to EeV

        #azimuth_sim = 180
        if primary!="proton":
            shower = {
        	   	"primary" : primary,        # primary (electron, pion)
        	   	"energy" : energy,               # EeV
        	   	"zenith" : zenith_sim,               # deg (GRAND frame)
        	   	"azimuth" : azimuth_sim,                # deg (GRAND frame)
        	   	"injection_height" : injection_height,    # m (injection height in the local coordinate system)
        	    "altitude" : GdAlt }   # m (alitude with respect to sealevel)

            # Perform the radiomorphing
            radiomorphing.process(sim_dir, shower, antennas, out_dir)

        datafile = file(inp_file)
        create_antpos(out_dir,datafile,slope)
    except:
        pass











########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
# Code line provided as an example by Anne Zilles.
'''
# Settings of the radiomorphing
data_dir = "/Users/nrenault/Desktop/GRAND/RadioMorphing/data/" #join(root_dir, "examples", "data")
sim_dir = join(data_dir, "GrandEventADetailed2") # folder containing your refernence shower simulations
out_dir = join(data_dir, "InterpolatedSignals2") # folder which will contain radio morphed traces afterwards
antennas = join(out_dir, "antpos.dat") # list of antenna positions you would like to simulate, stored in out_dir in the best case
height = 1500 #m above sea level

shower = {
    "primary" : "pion",        # primary (electron, pion)
    "energy" : 0.0691085,               # EeV
    "zenith" : 89.71,               # deg (GRAND frame)
    "azimuth" : 0.,                # deg (GRAND frame)
    "injection_height" : 1553.15,    # m (injection height in the local coordinate system)
    "altitude" : height }   # m (alitude with respect to sealevel)

# Perform the radiomorphing
radiomorphing.process(sim_dir, shower, antennas, out_dir)

shower = {
    "primary" : "pion",        # primary (electron, pion)
    "energy" : 0.644,               # EeV
    "zenith" : 89.3,               # deg (GRAND frame)
    "azimuth" : 0.,                # deg (GRAND frame)
    "injection_height" : 1846.44,    # m (injection height in the local coordinate system)
    "altitude" : GdAlt }   # m (alitude with respect to sealevel)


# Perform the radiomorphing
radiomorphing.process(sim_dir, shower, antennas, out_dir)
'''
