#!/usr/bin/env python

"""
                    Script PrepCR_zhaires_inp_files
                        Version 1.0
    Written by N. Renault-Tinacci from a script provided by C. Medina
                Using danton.py developped by V. Niess
"""

import json
import os, glob
import random
import shutil
import struct
import subprocess
import sys
import time
import numpy as np
import pylab as pl
import StringIO
import random
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
import modules

GEOMAGNET = (56.5, 63.18, 2.72) #Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg])
GdAlt=1500. #1500. #array altitude in m
Re= 6370949 # m, Earth radius

##########################################################################################################
def main():
    """ Main script allowing to produce the ZHAireS input files for each one of the showers from a DANTON library """

    ##################################################################
    # Test arguments
    if (len(sys.argv)<7 or len(sys.argv)>8):
        print """\
    This script will allow to produce the ZHAireS input files for cosmic-rays (hence down-going).
    It creates a regular rectangular array centered on the intersection point of the shower axis with the ground.
    The seed of each shower is uniformly randomized between 0 and 1.
    At the beginning of the script, you can set the altitude of the bottom of the array. By default GdAlt=1500 m
    The tables (pwr_EE, ZENITH and AZIMUTH) setting the energy, zenith and azimuth are, for now, hard coded and meant to be modified in the script.

    Inputs :
        work_dir = directory where all ZHAireS input files and the simulation results will be stored.
        slope = angle between horizontal and the slope along which the radio array is distributed [in degrees]
        height = maximum height that antennas distributed along the slope can reach [in m]
        step = separation between antennas [in m]
        Nx = number of antennas in the X direction
        Ny = number of antennas in the Y direction

    Ouput:
        The script will produce as many ZHAireS input files as there are combinations between the requested energies,
        zeniths and azimuths (respectively in tables pwr_EE, ZENITH and AZIMUTH).
        They will located in the work_dir+"/inp/" directory.

    Usage:  python PrepCR_zhaires.py work_dir slope height step Nx Ny

    Notes:
        The global variable DISPLAY allows you to turn on/off the display of the 3D map of the radio array and the print-out of the processed showers
        The "compute_antenna_pos" function can easily be modified/replaced to build a different antenna array configuration (star-shape pattern, circularly distributed, ...)
    """
    #in RASPASS North <=> azim=180 because azim and zenith reversed to trick Aires and then indicate dir of propag instead of dir of origin
        sys.exit(1)

    ##################################################################

    ##################################################################
    # Retrieve the input parameters
    ROOT = "/home/renault/CRs/"
    output_folder=str(sys.argv[1])
    DISPLAY = True #False #True

    ##################################################################
    # Array configuration parameters
    slope=float(sys.argv[2]) #slope [deg]
    hz=float(sys.argv[3]) #Array maximum height [m]
    sep=float(sys.argv[4]) #separation between antennas [m]
    nx = int(sys.argv[5]) #number of lines along x axis
    ny = int(sys.argv[6]) #number of lines along y axis

    ##################################################################
    #Showers parameters
    pwr_EE = np.array([18.],dtype=str) #np.array([17.5,18.,18.5,19.,19.5],dtype=str) #in eV
    AZIMUTH = np.array([0.],dtype=float) #np.array([0.,45.,90.,135.,180.],dtype=float) #in degrees in GRAND convention
    ZENITH = 180.-np.array([30.,65.,75.,85.],dtype=float) #180.-np.array([85.,80.,75.,70.,65.,60.],dtype=float) #in degrees in GRAND convention
    injh = 100e3 #above sea level in m

    ##################################################################
    print('******************************')
    print os.getcwd()

    Nit=0
    ANTENNAS=[]
    ANTENNAS=compute_antenna_pos(sep,nx,ny,GdAlt,slope)
    for azim in AZIMUTH:
        for zen in ZENITH:
            for eny in pwr_EE:
                task = 'EE1E'+eny+'-az'+str(int(azim))+'-zen'+str(int(zen))
                ANTENNAS2 = np.copy(ANTENNAS)

                a = np.sqrt((Re + injh)**2. - (Re+GdAlt)**2*np.sin(np.pi-np.deg2rad(zen))**2) - (Re+GdAlt)*np.cos(np.pi-np.deg2rad(zen))
                zen_inj = np.pi-np.arccos((a**2+(Re+injh)**2-Re**2)/(2*a*(Re+injh)))
                Xmax_primary,zen_inj2 = modules._getXmax('proton', 10**float(eny)/1e18, np.deg2rad(zen) # approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
                Xmax_height, Xmax_distance = modules._dist_decay_Xmax(zen_inj, injh, Xmax_primary) # d_prime: distance from decay point to Xmax

                ### Reduce the radio array to the shower geometrical footprint (we account for a footprint twice larger than the Cherenkov angle)
                CORE = random_array_pos(slope,sep) #
                ANTENNAS3 = reduce_antenna_array(Xmax_height,zen,azim,ANTENNAS2,CORE,DISPLAY)
                totito  = generate_input(task, 10**float(eny), azim, zen, ANTENNAS3)

                #Write file
                fileZha = output_folder+'inp/EE1E+'+eny+'-az'+str(int(azim))+'-zen'+str(int(zen))+'.inp'
                dir=os.path.dirname(fileZha)
                if not os.path.exists(dir):
                    os.makedirs(dir)
                inpfile = open(fileZha,"w+")
                inpfile.write(totito)
                inpfile.close()


##########################################################################################################
##########################################################################################################
###                                  Let's define useful functions                                     ###
##########################################################################################################
##########################################################################################################
def random_array_pos(slope=0.,sep=1e3):
    """ Compute a random offset for 2 or 3 of the space dimensions depending of the slope """

    CORE = np.array([0.,0.,0.])
    CORE[1] = random.uniform(-1.,1.)*sep/2 # always need an offset on Y (perpendicular to trajectory)
    if slope!=90. and slope!=0.: # random position along slope => x and z are random and their offset is related
        CORE[0] = random.uniform(0.,1.)*sep*np.cos(np.radians(slope))
        CORE[2] = CORE[0]*np.sin(np.radians(slope))
    elif slope==0.: # z = 0 => no offset in Z required
        CORE[0] = random.uniform(0.,1.)*sep
    elif slope==90.: # x = distance => no offset in X required
        CORE[2] = random.uniform(0,1.)*sep

    return CORE

##########################################################################################################
def getCerenkovAngle(h=100e3):
    """ Compute the Cherenkov angle of the shower at the altitude of injection """

    # h in meters
    n = 1.+325.e-6*np.exp(-0.1218*h*1e-3)      # Refractive index ZHAireS (see email M. Tueros 25/11/2016)
    alphac = np.arccos(1./n)
    return alphac

##########################################################################################################
def reduce_antenna_array(h=None,theta=None,phi=None,ANTENNAS=None,core=[0.,0.,0.],DISPLAY=False):
    """ Reduce the size of the initialized radio array to the shower geometrical footprint by computing the angle between shower and decay-point-to-antenna axes """
    """ theta = zenith in GRAND convention [in deg], h = Xmax height [in m] """

    zen,azim = GRANDtoZHAireS(theta,phi)
    zenr = np.radians(zen)
    azimr = np.radians(azim)
    ANTENNAS1 = np.copy(ANTENNAS)

    # Shift antenna array with the randomized core position
    ANTENNAS1[:,0] = ANTENNAS1[:,0]+core[0]
    ANTENNAS1[:,1] = ANTENNAS1[:,1]+core[1]
    ANTENNAS1[:,2] = ANTENNAS1[:,2]+core[2]

    # Compute angle between shower and decay-point-to-antenna axes
    #u_ant = ANTENNAS1-h*np.array([-np.tan(zenr),0.,1.],dtype=float)
    #u_sh = [np.sin(zenr),0.,-np.cos(zenr)]
    u_ant = ANTENNAS1-np.array([-(h-GdAlt)*np.tan(zenr)*np.cos(azimr),-(h-GdAlt)*np.tan(zenr)*np.sin(azimr),h],dtype=float)
    u_ant = (u_ant.T/np.linalg.norm(u_ant,axis=1)).T
    u_sh = np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), -np.cos(zenr)])
    ant_angle = np.arccos(np.matmul(u_ant, u_sh))

    # Remove antennas of the initial array that are located outside the "footprint"
    omegar = getCerenkovAngle(h)*5. #[in rad] # Accounting for a footprint four times larger than the Cherenkov angle
    angle_test = ant_angle<=omegar
    sel = np.where(angle_test)[0]
    ANTENNAS2 = ANTENNAS1[sel,:]

    # Remove the farthest antennas to reduce the number of antenna positions to simulate so that this number falls below 1000
    while np.shape(sel)[0]>999:
        #x_ant_max = np.max(np.abs(ANTENNAS2[:,0]))
        #antisel = np.where(np.abs(ANTENNAS2[:,0])==x_ant_max)[0]
        r_ant_max = np.max(np.sqrt(ANTENNAS2[:,0]**2+ANTENNAS2[:,1]**2))
        antisel = np.where(np.sqrt(ANTENNAS2[:,0]**2+ANTENNAS2[:,1]**2)==r_ant_max)[0]
        ANTENNAS2= np.delete(ANTENNAS2,antisel,0)
        sel= np.delete(sel,antisel,0)

    # 3D Display of the radio array
    if DISPLAY:
        ant_map_i = np.zeros(np.shape(ANTENNAS1)[0])
        ant_map_i[sel]=1.
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        array_display(ANTENNAS,ant_angle,'Shower axis to decay point-antenna axis angle map')
        array_display(ANTENNAS1,cc,'Selected antenna map')

    return ANTENNAS2

##########################################################################################################
def rotate_antenna_array(ANTENNAS=None,azim=0.):
    """ For a azimuth different of 0, one need to rotate the radio array with azimuth """
    """ so that the longest side of the array is aligned with shower axis """

    azimr = np.radians(azim)
    print azim
    if azimr>2*np.pi:
        azimr = azimr-2.*np.pi
    elif azimr<0.:
        azim = azim+2.*np.pi
    ANTENNAS2 = np.copy(ANTENNAS)
    ANTENNAS2[:,0] = ANTENNAS[:,0]*np.cos(azimr)-ANTENNAS[:,1]*np.sin(azimr)
    ANTENNAS2[:,1] = ANTENNAS[:,0]*np.sin(azimr)+ANTENNAS[:,1]*np.cos(azimr)
    ANTENNAS2[:,2] = ANTENNAS[:,2]

    # 3D Display of the radio array
    if DISPLAY:
        ant_map_i = np.ones(np.shape(ANTENNAS2)[0])
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        array_display(ANTENNAS2,cc,'Selected antenna map')

    return ANTENNAS2

##########################################################################################################
def compute_antenna_pos(step, nsidex,nsidey,GdAlt,inclin=0.):
    """ Generate antenna positions in a inclined plane @ a given distance from decay"""
    """ Return N positions (x,y,z) in Zhaires coordinates """
    xi,yi = (-0.5*nsidex*step+step*0.5,-0.5*nsidey*step+step*0.5)
    xf,yf = (0.5*nsidex*step+step*0.5,0.5*nsidey*step+step*0.5)
    xx, yy= np.meshgrid(np.arange(xi,xf,step), np.arange(yi,yf,step))
    zz=xx*np.tan(np.radians(inclin))
    xxr = np.ravel(xx)
    yyr = np.ravel(yy)
    zzr = np.ravel(zz)+GdAlt
    xyz = np.array([xxr, yyr, zzr]).T
    return xyz

##########################################################################################################
def generate_input(task, energy, azimuth, zenith, antennas):
    """Generate the input stream for ZHAIRES."""

    zen,azim = GRANDtoZHAireS(zenith,azimuth)

    seed = random.uniform(0.,1.)

    # Format the stream.
    stream = [

        "#########################",
        "TaskName {:s}".format(task),
        "PrimaryParticle proton",
        "PrimaryEnergy {:.2E} eV".format(energy),
        "PrimaryZenAngle {:.2f} deg".format(zenith),
        "PrimaryAzimAngle {:.2f} deg Magnetic".format(azimuth),
        "ForceModelName SIBYLL",
        #"SetGlobal RASPASSHeight {:.2f} m".format(injh),
        "RandomSeed {:.5f}".format(seed),
        "########################",
        "PropagatePrimary On",
        "SetGlobal RASPASSTimeShift 0.0",
        "SetGlobal RASPASSDistance 0.00"
    ]
    for a in antennas:
        stream.append("AddAntenna {:1.2f} {:1.2f} {:1.2f}".format(a[0],a[1],a[2]))

    stream += [
        "##########################",
        "TotalShowers 1",
        "RunsPerProcess Infinite",
        "ShowersPerRun 1",
        "Atmosphere 1",
        "AddSite Ulastai 42.55 deg 86.68 deg {:.3f} m".format(GdAlt),
        "Site Ulastai",
        "Date 1985 10 26",
        "GeomagneticField On",
        "GeomagneticField {:.4f} uT {:.2f} deg {:.2f} deg".format(*GEOMAGNET),
        "GroundAltitude {:.3f} m".format(GdAlt),
        "ObservingLevels 510 50 g/cm2   900 g/cm2",
        "PerShowerData Full",
        "SaveNotInFile lgtpcles All",
        "SaveNotInFile grdpcles All",
        "RLimsFile grdpcles 0.000001 m 10 km",
        "ResamplingRatio 100",
        "#########################",
        "RLimsTables 10 m 10 km",
        "ELimsTables 2 MeV 1 TeV",
        "ExportTables 5501 Opt a",
        "ExportTable 1205 Opt a",
        "ExportTable 1205 Opt as",
        "ExportTables 1293 Opt a",
        "ExportTables 1293 Opt as",
        "ExportTable 1793 Opt a",
        "ExportTable 1793 Opt as",
        "########################",
        "ForceLowEDecay Never",
        "ForceLowEAnnihilation Never",
        "########################",
        "ZHAireS On",
        "FresnelTime On",
        "FresnelFreq Off",
        "TimeDomainBin 1 ns",
        "AntennaTimeMin -100 ns",
        "AntennaTimeMax 1000 ns", #can be extended until 3e-6s but if then there is still nothing then there must be a problem somewhere
        "######################",
        "ElectronCutEnergy 1 MeV",
        "ElectronRoughCut 1 MeV",
        "GammaCutEnergy 1 MeV",
        "GammaRoughCut 1 MeV",
        "ThinningEnergy 1.e-4 Relative", #It can be 1e-5, 1e-6 or below. But running time inversely proportional to it.
        "ThinningWFactor 0.06"
    ]
    return "\n".join(stream)

##########################################################################################################
def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)

##########################################################################################################
def array_display(ANTENNAS=None,datamap=None,title=None):
    if len(ANTENNAS[:,0])!=0:
        fig1 = pl.figure(1,figsize=(5*3.13,3.8*3.13))
        binmap = ListedColormap(['white', 'black'], 'indexed')
        dar=(np.max(ANTENNAS[:,0])-np.min(ANTENNAS[:,0]))/(np.max(ANTENNAS[:,1])-np.min(ANTENNAS[:,1]))
        if dar==0:
            dar=1
        xlbl='X [m]'
        ylbl='Y [m]'
        zlbl='Z [m]'

        ax = pl.gca(projection='3d')
        ax.scatter(ANTENNAS[:,0]*1.,ANTENNAS[:,1],ANTENNAS[:,2],c=datamap)
        ax.set_title(title)
        ax.view_init(25,-130)
        pl.xlabel(xlbl)
        pl.ylabel(ylbl)
        ax.set_zlabel(zlbl)
        #pl.gca().set_aspect(1,adjustable='box')
        pl.gca().set_aspect('equal')

        pl.show()
    return

##########################################################################################################
def GRANDtoZHAireS(zen_DANTON=None, azim_DANTON=0):
    """ Convert coordinates from DANTON convention to ZHAireS convention """

    zen = 180. - zen_DANTON
    azim = azim_DANTON - 180.
    if azim>360:
        azim = azim-360.
    elif azim<0.:
        azim = azim+360.
    return [zen,azim]

##########################################################################################################
def ZHAireStoGRAND(zen_ZHAireS=None, azim_ZHAireS=0):
    """ Convert coordinates from DANTON convention to ZHAireS convention """

    zen = (180.-zen_ZHAireS)
    azim = azim_ZHAireS + 180.
    if azim>360:
        azim = azim-360.
    elif azim<0.:
        azim = azim+360.
    return [zen,azim]

##########################################################################################################
def DANTONtoGRAND(zen_DANTON=None, azim_DANTON=0):
    """ Convert coordinates from DANTON convention to ZHAireS convention """

    zen = zen_DANTON
    azim = azim_DANTON
    if azim>360:
        azim = azim-360.
    elif azim<0.:
        azim = azim+360.
    return zen,azim

##########################################################################################################
if __name__ == '__main__':
    main()
