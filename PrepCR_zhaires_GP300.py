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
import matplotlib #.pyplot as plt
matplotlib.use('Agg')
import pylab as pl
import StringIO
import random
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('/home/renault/scripts_clean/')
import modules
sys.path.append('/home/renault/retro/lib/')
sys.path.append('/home/renault/retro/lib/python/')
from grand_tour import Topography
sys.path.append('/home/renault/retro/deps/turtle/src')
import turtle
import os
os.chdir('/home/renault/retro/')

GEOMAGNET = (56.5, 63.18, 2.72) #Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg])
Re= 6370949 # m, Earth radius
MIN_ANTENNAS = 8
latitude, longitude = 43.5,94.0
#topo_path = 'share/topography'
topo_path = '/data75/renault/topography'
topo = Topography(latitude=latitude, longitude=longitude,path=topo_path, stack_size=121)
GdAlt = 2016.41083928 #m = alt of the lowest antenna of GP300 #topo.ground_altitude(0., 0.)

##########################################################################################################
def main():
    """ Main script allowing to produce the ZHAireS input files for each one of the showers from a DANTON library """

    ##################################################################
    # Test arguments
    if (len(sys.argv)<2 or len(sys.argv)>3):
        print """\
    This script will allow to produce the ZHAireS input files for cosmic-rays (hence down-going).
    It reads the antenna positions defined for GRANDproto300 (x,y,altitude). The energy is hardcoded.
    AZIMUTH [0,360[, ZENITH [90, 140[, core position (chosen between -draw_area_size and +draw_area_size around the origin both in x and y directions) are randomized.
    Then a selection of the antenna falling into the Cherenkov ring is done and the process is iterated until  more than MIN_ANTENNAS antennas are contained in the Cherenkov ring. The number of tries needed to meet
    the criterion to get one event is recorded together with a bunch of informations about the shower and the array.

    Inputs :
        work_dir = directory where all ZHAireS input files and the simulation results will be stored.

    Ouput:
        The script will produce as many ZHAireS input files as there are combinations between the requested energies,
        zeniths and azimuths (respectively in tables pwr_EE, ZENITH and AZIMUTH).
        They will located in the work_dir+"/inp/" directory.

    Usage:  python PrepCR_zhaires_GP300.py work_dir

    Notes:
        The global variable DISPLAY allows you to turn on/off the display of the 3D map of the radio array and the print-out of the processed showers
        The "compute_antenna_pos" function can easily be modified/replaced to build a different antenna array configuration (star-shape pattern, circularly distributed, ...)
    """
    #in RASPASS North <=> azim=180 because azim and zenith reversed to trick Aires and then indicate dir of propag instead of dir of origin
        sys.exit(1)
    print('******************************')

    ##################################################################
    # Retrieve the input parameters
#    ROOT = "/home/renault/CRs_GP300/"
    output_folder=str(sys.argv[1])
    DISPLAY = False

    ##################################################################
    #Showers parameters
    pwr_EE = np.array([16.5,17.,17.5,18.,18.5,19.,19.5],dtype=str) #in eV
    injh = 100e3 #above sea level in m
    draw_area_size = 50e3 #m

    ##################################################################
    # Compute an antenna array
    #ANTENNAS = np.loadtxt('/Users/nrenault/Downloads/GP300/5_good/GP300_antpos.txt')
    #antcoord = np.loadtxt('/Users/nrenault/Downloads/GP300/5_good/GP300_antcoord.txt')
    antcoord = np.loadtxt('/home/renault/CRs_GP300/GP300_layout/5_good/GP300_antcoord.txt')
    ANTENNAS = np.loadtxt('/home/renault/CRs_GP300/GP300_layout/5_good/GP300_antpos.txt')
    ANTENNAS[:,2] = antcoord[:,2]

    ### Reduce the radio array to the shower geometrical footprint (we account for a footprint twice larger than the Cherenkov angle)
    dx = draw_area_size #+ (np.amax(ANTENNAS[:,0])-np.amin(ANTENNAS[:,0]))/2
    dy = draw_area_size #+ (np.amax(ANTENNAS[:,1])-np.amin(ANTENNAS[:,1]))/2

    shower_infos = []
    Nwish = 500
    ievt = 0
    for eny in pwr_EE: #[:1]:
        Nevents = 0
        while Nevents<Nwish:
            ievt +=1
            azim = random.uniform(0.,1.)*360. #GRAND convention
            zen = random.uniform(0.,1.)*50.+90. #GRAND convention #zenith between 90 and 140 deg

            # Compute informations on Xmax
            Xmax_primary,zen_inj = modules._getXmax('proton', 10**float(eny)/1e9, np.deg2rad(zen),GdAlt,injh)
            # approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
            Xmax_height, Xmax_distance = modules._dist_decay_Xmax(zen_inj, injh, Xmax_primary)
            # d_prime: distance from decay point to Xmax

            # Sample the shower core position.
            b = np.mean(ANTENNAS[:,0:], axis=0)
            ANTENNAS1 = np.copy(ANTENNAS)
            x = random.uniform(b[0] - dx, b[0] + dx)
            y = random.uniform(b[1] - dy, b[1] + dy)

            try:
                altcore = topo.ground_altitude(x, y)
            except:
                altcore = 2000.

            ANTENNAS2 = reduce_antenna_array(Xmax_height,zen,azim,ANTENNAS1,x,y,altcore,10**float(eny)/1e9,DISPLAY) #Xmax_height
            if len(ANTENNAS2) >= MIN_ANTENNAS:
                Nevents += 1
                task = 'EE1E'+eny+'_evt'+str(Nevents)
                print Nevents,ievt,azim,zen,x,y,len(ANTENNAS2)
                totito  = generate_input(task, 10**float(eny), azim, zen, ANTENNAS2-[x,y,0.])

                #Write file
                fileZha = output_folder+'inp/'+task+'.inp'
                dir=os.path.dirname(fileZha)
                if not os.path.exists(dir):
                    os.makedirs(dir)
                    os.makedirs(dir+'/fig/')
                inpfile = open(fileZha,"w+")
                inpfile.write(totito)
                inpfile.close()

                plot_drawn_pos(ANTENNAS1,ANTENNAS2,[x,y,altcore],output_folder,task,DISPLAY)
                shower_infos.append([ievt,Nevents,10**float(eny)/1e18,azim,zen,len(ANTENNAS2),b[0],b[1],x,y,altcore,dx,dy,GdAlt,MIN_ANTENNAS])
            else:
                shower_infos.append([ievt,0,10**float(eny)/1e18,azim,zen,len(ANTENNAS2),b[0],b[1],x,y,altcore,dx,dy,GdAlt,MIN_ANTENNAS])

    shower_info_file = output_folder+'/inp/summary_table.txt'
    shower_infos = np.array(shower_infos)
    hdr='\n'.join(["NumEvt Nevents_selected[0 if not selected] Energy[EeV] Azimuth_GRAND[deg] Zenith_GRAND[deg] NAntennas Xmean_array[m] Ymean_array[m] Xcore[m] Ycore[m] AltCore[m] dx[m] dy[m] GroundAltitudeAboveSeaLevel[m] MinAntennas"])
    np.savetxt(shower_info_file,shower_infos,fmt='%d   %d   %0.2f   %0.2f   %0.2f   %d   %0.2f   %0.2f   %0.2f   %0.2f   %0.2f   %0.1f   %0.1f %0.2f %d' ,header=hdr)

##########################################################################################################
##########################################################################################################
###                                  Let's define useful functions                                     ###
##########################################################################################################
##########################################################################################################
def plot_drawn_pos(ANTENNAS_i,ANTENNAS_sel,pos_core,output_folder,task,DISPLAY):
    #ANTENNAS_i are the antenna positions corrected by the position of the intersection between the shower axis and the ground
    #ANTENNAS_sel are the selected antenna positions
    #CORE is the random position

    xlbl='X [m]'
    ylbl='Y [m]'
    zlbl='Z [m]'

    #2D (x-y)
    fig, ax = pl.subplots()
    pl.scatter(ANTENNAS_i[:,0],ANTENNAS_i[:,1],c='b',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(ANTENNAS_sel[:,0],ANTENNAS_sel[:,1],c='g',edgecolors='none') #arrays are move back around (0,0,GdAlt)
    pl.scatter(pos_core[0],pos_core[1],c='r',edgecolors='none')
    pl.scatter(0.,0.,c='y',edgecolors='none')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    figname = output_folder+'/inp/fig/plots_'+task+'.png'
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
def reduce_antenna_array(h,theta,phi,ANTENNAS,x,y,altcore,Esh,DISPLAY=False):
    """ Reduce the size of the initialized radio array to the shower geometrical footprint by computing the angle between shower and decay-point-to-antenna axes """
    """ theta = zenith in GRAND convention [in deg], h = Xmax height [in m] """

    zenr = np.radians(theta)
    azimr = np.radians(phi)
    ANTENNAS1 = np.copy(ANTENNAS)

    # Compute angle between shower and decay-point-to-antenna axes
    Xmax_pos = np.array([x+(h-altcore)*np.tan(zenr)*np.cos(azimr),y+(h-altcore)*np.tan(zenr)*np.sin(azimr),h],dtype=float)
    u_sh = np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), np.cos(zenr)])

    # Remove antennas of the initial array that are located outside the "footprint"
    sel = select(ANTENNAS1,Esh,Xmax_pos, u_sh,topo)
    ANTENNAS2 = ANTENNAS1[sel,:]

    # 3D Display of the radio array
    if DISPLAY and len(ANTENNAS2)>MIN_ANTENNAS:
        ant_map_i = np.zeros(np.shape(ANTENNAS1)[0])
        ant_map_i[sel]=1.
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        #array_display(ANTENNAS,ant_angle,'Shower axis to decay point-antenna axis angle map')
        array_display(ANTENNAS1,cc,'Selected antenna map')
    return ANTENNAS2

##########################################################################################################
def select(ra, shower_energy, position, direction,topography,check_xmax=True, shadowing=True):
    deltar = 200.

    # Cone parameters
    cherenkov_angle = 3. #1.5
    print 'Cherenkov angle = ',cherenkov_angle
    gamma = np.deg2rad(cherenkov_angle)
    zcmin = 14E+03  # m
    zcmax = 165E+03 * shower_energy / 1E+09 + 55E+03  # m

    # Check if the shower crashes into a mountain early, before xmax
    if check_xmax:
        s = np.arange(0., zcmin + deltar, deltar)
        xs, ys, zs = [position[i] + direction[i] * s for i in xrange(3)]
        zg = [topography.ground_altitude(xi, yi) for xi, yi in zip(xs, ys)]
        if (zs <= zg).any():
            return []

    # Select the antenna(s) within the cone
    dr = ra[:, :3] - position
    zp = np.dot(dr, direction)
    rp2 = np.sum(dr**2, axis=1) - zp**2
    test_radius = rp2 <= ((zp - zcmin) * np.tan(gamma))**2
    test_edgemin = zp >= zcmin
    test_edgemax = zp <= zcmax
    index = np.nonzero(test_radius & test_edgemin & test_edgemax)[0]
    if len(index) == 0:
        return index

    if not shadowing:
        return ra[index, :].tolist()

    # Check for shadowing
    r0 = position + zcmin * np.array(direction)

    def check_shadowing(i):
        u = ra[i, :3] - r0
        d = np.linalg.norm(u)
        u /= d
        s = np.arange(0., d, deltar)
        xj, yj, zj = [r0[j] + u[j] * s for j in xrange(3)]
        zg = [topography.ground_altitude(x, y) for x, y in zip(xj, yj)]
        if (zj <= zg).any():
            return False
        return True

    K = filter(check_shadowing, index)
    return K #ra[K, :].tolist(),K

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
        "AddSite Bailikun 43.61 deg 93.82 deg {:.3f} m".format(GdAlt),
        "Site Bailikun",
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
        ax.view_init(elev=90., azim=-90.)
        #pl.gca().set_aspect(1,adjustable='box')
        #pl.gca().set_aspect('equal')
        #pl.axis('equal')

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
