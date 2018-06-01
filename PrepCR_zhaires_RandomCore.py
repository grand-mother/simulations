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
#import matplotlib #.pyplot as plt
#matplotlib.rcParams['backend'] = "WXAgg"
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as pl
import pylab as pl
import StringIO
import random
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
import modules

GEOMAGNET = (56.5, 63.18, 2.72) #Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg])
GdAlt=1500. #array altitude in m
Re= 6370949 # m, Earth radius
MIN_ANTENNAS = 5

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
    If the direction and energy are set, the intersection point between the shower axis and ground is randomly chosen between
    -draw_area_size and +draw_area_size around the array edges both in x and y directions. Then a selection of the antenna falling into the Cherenkov ring is done,
    and the process is iterated until  more than MIN_ANTENNAS antennas are contained in the Cherenkov ring. The number of tries needed to meet
    the criterion to get one event is recorded together with a bunch of informations about the shower and the array.

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

    Usage:  python PrepCR_zhaires_RandomCore.py work_dir slope height step Nx Ny

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
    DISPLAY = False #True #

    ##################################################################
    # Array configuration parameters
    slope=float(sys.argv[2]) #slope [deg]
    hz=float(sys.argv[3]) #Array maximum height [m]
    sep=float(sys.argv[4]) #separation between antennas [m]
    nx = int(sys.argv[5]) #number of lines along x axis
    ny = int(sys.argv[6]) #number of lines along y axis

    ##################################################################
    #Showers parameters
    pwr_EE = np.array([17.5,18.,18.5,19.,19.5],dtype=str) #in eV
    AZIMUTH = np.array([0.,45.,90.,135.,180.],dtype=float) #in degrees in GRAND convention
    ZENITH = 180.-np.array([85.,80.,75.,70.,65.,60.],dtype=float) #in degrees in GRAND convention
    injh = 100e3 #above sea level in m
    draw_area_size = 50e3

    ##################################################################
    print('******************************')
    print os.getcwd()

    # Compute an antenna array
    ANTENNAS=[]
    ANTENNAS=compute_antenna_pos(sep,nx,ny,GdAlt,slope)
    print ANTENNAS
    print np.amax(ANTENNAS[:,0]),np.amin(ANTENNAS[:,0])
    print np.amax(ANTENNAS[:,1]),np.amin(ANTENNAS[:,1])
    exit()
    ### Reduce the radio array to the shower geometrical footprint (we account for a footprint twice larger than the Cherenkov angle)
    CORE = random_array_pos(slope,sep) #[0.,0.,0.] #
    ANTENNAS[:,0] = ANTENNAS[:,0]+CORE[0]
    ANTENNAS[:,1] = ANTENNAS[:,1]+CORE[1]
    ANTENNAS[:,2] = ANTENNAS[:,2]+CORE[2]
    dx = draw_area_size + (np.amax(ANTENNAS[:,0])-np.amin(ANTENNAS[:,0]))/2
    dy = draw_area_size + (np.amax(ANTENNAS[:,1])-np.amin(ANTENNAS[:,1]))/2

    shower_infos = []
    Nwish = 20
    for ievt in range(0,Nwish):
        for azim in AZIMUTH:
            for zen in ZENITH:
                for eny in pwr_EE:
                    task = 'EE1E'+eny+'-az'+str(int(azim))+'-zen'+str(int(zen))+'_evt'+str(ievt)

                    # Compute informations on Xmax
                    Xmax_primary,zen_inj = modules._getXmax('proton', 10**float(eny)/1e18, np.deg2rad(zen)) # approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
                    Xmax_height, Xmax_distance = modules._dist_decay_Xmax(zen_inj, injh, Xmax_primary) # d_prime: distance from decay point to Xmax

                    # Sample the shower core position.
                    b = np.mean(ANTENNAS[:,0:], axis=0)

                    pos_tab= []
                    Nevents = 0
                    while True:
                        Nevents += 1
                        ANTENNAS2 = np.copy(ANTENNAS)
                        x = random.uniform(b[0] - dx, b[0] + dx)
                        y = random.uniform(b[1] - dy, b[1] + dy)
                        ANTENNAS2 = ANTENNAS2-[x,y,0]
                        pos_tab.append([x,y])
                        if Nevents %1e5==0:
                            random.seed()
                            print Nevents
                        ANTENNAS3 = reduce_antenna_array(Xmax_height,zen,azim,ANTENNAS2,DISPLAY)
                        if len(ANTENNAS3) >= MIN_ANTENNAS:
                            break
                    pos_tab = np.array(pos_tab)
                    totito  = generate_input(task, 10**float(eny), azim, zen, ANTENNAS3)

                    #Write file
                    fileZha = output_folder+'inp/'+task+'.inp'
                    dir=os.path.dirname(fileZha)
                    if not os.path.exists(dir):
                        os.makedirs(dir)
                        os.makedirs(dir+'/fig/')
                    inpfile = open(fileZha,"w+")
                    inpfile.write(totito)
                    inpfile.close()

                    plot_drawn_pos(ANTENNAS2,ANTENNAS3,pos_tab,output_folder,task,DISPLAY)
                    print 'Nevents =',Nevents
                    shower_infos.append([ievt,10**float(eny)/1e18,azim,zen,Nevents,len(ANTENNAS3),x,y,dx,dy,sep,nx,ny,slope,hz,GdAlt,MIN_ANTENNAS])
    shower_info_file = output_folder+'/inp/summary_table.txt'
    shower_infos = np.array(shower_infos)
    hdr='\n'.join(["NumEvt Energy[EeV] Azimuth_GRAND[deg] Zenith_GRAND[deg] Ntry NAntennas x[m] y[m] dx[m] dy[m] step[m] Nx Ny Slope[deg] MountainHeight[m] GroundAltitudeAboveSeaLevel[m] MinAntennas", "xmin= "+str(np.amin(ANTENNAS[:,0]))+" xmax= "+str(np.amax(ANTENNAS[:,0]))+" ymin= "+str(np.amin(ANTENNAS[:,1]))+" ymax= "+str(np.amax(ANTENNAS[:,1]))])
    np.savetxt(shower_info_file,shower_infos,fmt='%d %0.2f   %0.1f   %0.1f   %d   %d   %0.2f   %0.2f   %0.2f   %0.2f   %0.1f   %d   %d   %0.1f   %0.1f   %0.1f   %d' ,header=hdr)

##########################################################################################################
##########################################################################################################
###                                  Let's define useful functions                                     ###
##########################################################################################################
##########################################################################################################
def plot_drawn_pos(ANTENNAS_i,ANTENNAS_sel,pos_tab,output_folder,task,DISPLAY):
    #ANTENNAS_i are the antenna positions corrected by the position of the intersection between the shower axis and the ground
    #ANTENNAS_sel are the selected antenna positions
    #CORE is the random position
    fig, ax = pl.subplots()
    pl.scatter(ANTENNAS_i[:,0]+pos_tab[-1,0],ANTENNAS_i[:,1]+pos_tab[-1,1],c='b',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(ANTENNAS_sel[:,0]+pos_tab[-1,0],ANTENNAS_sel[:,1]+pos_tab[-1,1],c='g',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(pos_tab[:,0],pos_tab[:,1],c='k',edgecolors='none')
    pl.scatter(pos_tab[-1,0],pos_tab[-1,1],c='r',edgecolors='none')
    figname = output_folder+'/inp/fig/plots_'+task+'.png'
    pl.savefig(figname,dpi=350)
    if DISPLAY:
        pl.show()
    pl.close()

    fig, ax = pl.subplots()
    pl.scatter(ANTENNAS_i[:,0],ANTENNAS_i[:,1],c='b',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(ANTENNAS_sel[:,0],ANTENNAS_sel[:,1],c='g',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(0.,0.,c='r',edgecolors='none')
    figname = output_folder+'/inp/fig/plots_'+task+'_real.png'
    pl.savefig(figname,dpi=350)
    if DISPLAY:
        pl.show()
    pl.close()
    return

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
def reduce_antenna_array(h=None,theta=None,phi=None,ANTENNAS=None,DISPLAY=False):
    """ Reduce the size of the initialized radio array to the shower geometrical footprint by computing the angle between shower and decay-point-to-antenna axes """
    """ theta = zenith in GRAND convention [in deg], h = Xmax height [in m] """

    zenr = np.radians(theta)
    azimr = np.radians(phi)
    ANTENNAS1 = np.copy(ANTENNAS)

    # Compute angle between shower and decay-point-to-antenna axes
    Xmax_pos = np.array([(h-GdAlt)*np.tan(zenr)*np.cos(azimr),(h-GdAlt)*np.tan(zenr)*np.sin(azimr),h],dtype=float)
    u_ant = ANTENNAS1-Xmax_pos
    u_ant = (u_ant.T/np.linalg.norm(u_ant,axis=1)).T
    u_sh = np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), np.cos(zenr)])
    ant_angle = np.arccos(np.matmul(u_ant, u_sh))

    # Remove antennas of the initial array that are located outside the "footprint"
    omegar = getCerenkovAngle(h)*4. #[in rad] # Accounting for a footprint four times larger than the Cherenkov angle
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
    if DISPLAY and len(sel)>MIN_ANTENNAS:
        ant_map_i = np.zeros(np.shape(ANTENNAS1)[0])
        ant_map_i[sel]=1.
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        #array_display(ANTENNAS,ant_angle,'Shower axis to decay point-antenna axis angle map')
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
        "PrimaryZenAngle {:.2f} deg".format(zen),
        "PrimaryAzimAngle {:.2f} deg Magnetic".format(azim),
        "ForceModelName SIBYLL",
        #"SetGlobal RASPASSHeight {:.2f} m".format(injh),
        "RandomSeed {:.5f}".format(seed),
        "########################",
        "PropagatePrimary On",
        "SetGlobal RASPASSTimeShift 0.0",
        "SetGlobal RASPASSDistance 0.00"
    ]
    for a in antennas:
        a[2] = a[2]-GdAlt
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
