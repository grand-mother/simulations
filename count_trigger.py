#!/usr/bin/env python
#import matplotlib
#matplotlib.use('Agg')

import os, glob
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys
import numpy as np
import linecache
from matplotlib import rc
import scipy.interpolate as itp
from matplotlib.image import NonUniformImage

ROOT = "/Users/nrenault/Desktop/GRAND/ToyModel/"

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<1 or len(sys.argv)>1):
    print """\
        This script will compute and plot the trigger efficiency for a set of simulation (designed at the moment for the ToyModel)

        Usage:  python count_trigger.py

        Note: all parameters are hard coded at the moment.
        """
    sys.exit(1)
###########################

Dd=np.array([20000,30000,40000,60000,80000,100000]) #distance to decay point in m
Dd2=[15000,25000,35000,50000,70000,90000,100000]
alfa=np.array([0,5,10,15,20,45,90]) #slope angle in degrees
alfa2=[-2.5,2.5,7.5,12.5,17.5,32.5,67.5,90.]

hz=3000 #Mountain Height
threshold_cons=150. #noise=15muV => case with 10*noise
threshold_aggr=50. #noise=15muV => case with 10/3*noise

Nserie=2
sep=[500,500,500,250,250,250,250] #separation between antennas
Nant_threshold=[32,32,32,128,128,128,128] #related to the step (=sep)
EE='5e+17' #1E+10' #Neutrino energy e
flag_val = False #True # If we want to plot valentin results or not
discarded_showers = np.array([],dtype=str) #1002503,1005547,100754,1011947,1015623,1029503,1030112,1031492,1032898,1036059,1041419,1042706,1045873,1047816,1082567],dtype=str)
#suffix = '_50-200MHz'
suffix = ''

####################################################################################
Nconfig_ew_cons = np.zeros((np.size(Dd),np.size(alfa)),dtype=float)
Nconfig_ns_cons = np.copy(Nconfig_ew_cons)
Nconfig_tot_cons = np.copy(Nconfig_ew_cons)
Nconfig_ew_aggr = np.copy(Nconfig_ew_cons)
Nconfig_ns_aggr = np.copy(Nconfig_ew_cons)
Nconfig_tot_aggr = np.copy(Nconfig_ew_cons)
for iDd in range(np.size(Dd)):
    for ialfa in range(np.size(alfa)):
        if alfa[ialfa]>0 or (alfa[ialfa]==0 and Dd[iDd]==20000):

            dist = []
            alpha = []
            showerID = []
            New_cons = []
            Nns_cons = []
            Nup_cons = []
            Ntot_cons = []
            New_aggr = []
            Nns_aggr = []
            Nup_aggr = []
            Ntot_aggr = []
            Vp2p_ew = []
            Vp2p_ns = []
            Vp2p_up = []
            Vp2p_tot = []
            for iserie in range(1,Nserie+1):
                test_file = ROOT+'detection/detection_count_EE'+EE+'_D'+str(Dd[iDd])+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+suffix+'.txt'
                dist_tmp, alpha_tmp, showerID_tmp, New_cons_tmp, Nns_cons_tmp, Nup_cons_tmp, Ntot_cons_tmp, New_aggr_tmp, Nns_aggr_tmp, Nup_aggr_tmp, Ntot_aggr_tmp, Vp2p_ew_tmp, Vp2p_ns_tmp, Vp2p_up_tmp, Vp2p_tot_tmp = np.loadtxt(test_file,comments='#',dtype=str,usecols=range(15),unpack=True,skiprows=1)
                dist = np.concatenate((dist,dist_tmp))
                alpha = np.concatenate((alpha,alpha_tmp))
                showerID = np.concatenate((showerID,showerID_tmp))
                New_cons = np.concatenate((New_cons,New_cons_tmp))
                Nns_cons = np.concatenate((Nns_cons,Nns_cons_tmp))
                Nup_cons = np.concatenate((Nup_cons,Nup_cons_tmp))
                Ntot_cons = np.concatenate((Ntot_cons,Ntot_cons_tmp))
                New_aggr = np.concatenate((New_aggr,New_aggr_tmp))
                Nns_aggr = np.concatenate((Nns_aggr,Nns_aggr_tmp))
                Nup_aggr = np.concatenate((Nup_aggr,Nup_aggr_tmp))
                Ntot_aggr = np.concatenate((Ntot_aggr,Ntot_aggr_tmp))

                Vp2p_ew = np.concatenate((Vp2p_ew,Vp2p_ew_tmp))
                Vp2p_ns = np.concatenate((Vp2p_ns,Vp2p_ns_tmp))
                Vp2p_up = np.concatenate((Vp2p_up,Vp2p_up_tmp))
                Vp2p_tot = np.concatenate((Vp2p_tot,Vp2p_tot_tmp))

            showerID = np.array(showerID)
            New_cons=np.array(New_cons,dtype=float)
            Nns_cons=np.array(Nns_cons,dtype=float)
            Ntot_cons=np.array(Ntot_cons,dtype=float)
            New_aggr=np.array(New_aggr,dtype=float)
            Nns_aggr=np.array(Nns_aggr,dtype=float)
            Ntot_aggr=np.array(Ntot_aggr,dtype=float)

            if np.size(discarded_showers)>0:
                ind_ID = np.zeros((len(discarded_showers),1))
                for ii in range(len(discarded_showers)):
                    ind_ID[ii] = np.where(showerID==discarded_showers[ii])
                New_cons=np.delete(New_cons,ind_ID)
                Nns_cons=np.delete(Nns_cons,ind_ID)
                Ntot_cons=np.delete(Ntot_cons,ind_ID)
                New_aggr=np.delete(New_aggr,ind_ID)
                Nns_aggr=np.delete(Nns_aggr,ind_ID)
                Ntot_aggr=np.delete(Ntot_aggr,ind_ID)

            Nconfig_ew_cons[iDd,ialfa] = np.size(np.where(New_cons>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))
            Nconfig_ns_cons[iDd,ialfa] = np.size(np.where(Nns_cons>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))
            Nconfig_tot_cons[iDd,ialfa] = np.size(np.where(Ntot_cons>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))
            Nconfig_ew_aggr[iDd,ialfa] = np.size(np.where(New_aggr>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))
            Nconfig_ns_aggr[iDd,ialfa] = np.size(np.where(Nns_aggr>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))
            Nconfig_tot_aggr[iDd,ialfa] = np.size(np.where(Ntot_aggr>=Nant_threshold[ialfa]))/float(np.size(Ntot_aggr))

####################################################################################
#Read and plot the results for Valentin = cone model parametrization
if flag_val:
    valentin_aggr_file = '/Users/nrenault/Desktop/GRAND/ToyModel/agressive_inter_'+EE+'eV.txt'
    valentin_cons_file = '/Users/nrenault/Desktop/GRAND/ToyModel/conservative_inter_'+EE+'eV.txt'
    data_geom_cons = np.loadtxt(valentin_cons_file)
    data_geom_aggr = np.loadtxt(valentin_aggr_file)
    detection_rate_cons = []
    detection_rate_aggr = []
    antenna_criteria = 128
    slope = np.unique(data_geom_cons[:,1])
    distance = np.unique(data_geom_cons[:,0])
    for d in distance:
        a_cons = data_geom_cons[np.logical_or.reduce([data_geom_cons[:,0]==d])]
        a_aggr = data_geom_aggr[np.logical_or.reduce([data_geom_aggr[:,0]==d])]
        for s in slope:
            b_cons = a_cons[np.logical_or.reduce([a_cons[:,1]==s])]
            b_aggr = a_aggr[np.logical_or.reduce([a_aggr[:,1]==s])]
            detection_cons = np.sum(b_cons[:,4]>=antenna_criteria)
            detection_aggr = np.sum(b_aggr[:,4]>=antenna_criteria)
            detection_rate_cons.append([d,s,detection_cons/float(len(b_cons[:,4]))])
            detection_rate_aggr.append([d,s,detection_aggr/float(len(b_aggr[:,4]))])

    detection_rate_aggr = np.array(detection_rate_aggr)
    detection_rate_cons = np.array(detection_rate_cons)

####################################################################################
##### Plot #####
#pl.rc('text', usetex=True)
#pl.rc('font', family='serif')
DISPLAY_Dd=1
DISPLAY_alfa=1
DISPLAY_2D=1
save_folder=ROOT+'detection/'
fig1 = pl.subplots(nrows=2, ncols=2,figsize=(4.*3.13,3.1*3.13))

#Plot vs Dd
if DISPLAY_Dd==1:
    pl.subplot(321)
    ind=np.where(alfa==10)[0][0]

    pl.plot(Dd/1e3,Nconfig_ew_cons[:,ind],label='Conservative (Vew)',linewidth=3,linestyle='--')
    #pl.plot(Dd/1e3,Nconfig_ns_cons[:,ind],label='Conservative (Vns)',linewidth=3,linestyle='--')
    #pl.plot(Dd/1e3,Nconfig_tot_cons[:,ind],label='Conservative (Vtot)',linewidth=3,linestyle='--')

    pl.plot(Dd/1e3,Nconfig_ew_aggr[:,ind],label='Aggressive (Vew)',linewidth=3,linestyle='-')
    #pl.plot(Dd/1e3,Nconfig_ns_aggr[:,ind],label='Aggressive (Vns)',linewidth=3,linestyle='-')
    #pl.plot(Dd/1e3,Nconfig_tot_aggr[:,ind],label='Aggressive (Vtot)',linewidth=3,linestyle='-')

    if flag_val:
        ind_val = np.where(detection_rate_aggr[:,1]==10)
        aaa_aggr = detection_rate_aggr[ind_val,0][0]/1e3
        bbb_aggr = detection_rate_aggr[ind_val,2][0]
        aaa_cons = detection_rate_cons[ind_val,0][0]/1e3
        bbb_cons = detection_rate_cons[ind_val,2][0]
        pl.plot(aaa_cons,bbb_cons,label='Conservative (Cone)',linewidth=3,linestyle='--')
        pl.plot(aaa_aggr,bbb_aggr,label='Aggressive (Cone)',linewidth=3,linestyle='-')

    pl.title('alpha = '+str(alfa[ind])+' deg')
    pl.xlabel('Distance (km)')
    #pl.ylabel(r'$\frac{Ntrig(<128antennas)}{Nsimulated showers(135)}$',)
    #pl.ylabel('Ntrig(<128antennas) / Nsimu showers(=135)',)
    pl.ylabel('Trigger rate')
    pl.legend(loc='best',fontsize=10)
    #pl.legend(bbox_to_anchor=(0., -0.3, 0.75, .102), loc=2, ncol=2, mode="expand", borderaxespad=0.,fontsize=10)

#Plot vs alpha
if DISPLAY_alfa==1:
    pl.subplot(322)
    ind2=np.where(Dd==30000)[0][0]

    pl.plot(alfa,Nconfig_ew_cons[ind2,:],label='Conservative (Vew)',linewidth=2,linestyle='--')
    #pl.plot(alfa,Nconfig_ns_cons[ind2,:],label='Conservative (Vns)',linewidth=2,linestyle='--')
    #pl.plot(alfa,Nconfig_tot_cons[ind2,:],label='Conservative (Vtot)',linewidth=2,linestyle='--')

    pl.plot(alfa,Nconfig_ew_aggr[ind2,:],label='Aggressive (Vew)',linewidth=2,linestyle='-')
    #pl.plot(alfa,Nconfig_ns_aggr[ind2,:],label='Aggressive (Vns)',linewidth=2,linestyle='-')
    #pl.plot(alfa,Nconfig_tot_aggr[ind2,:],label='Aggressive (Vtot)',linewidth=2,linestyle='-')

    if flag_val:
        ind_val = np.where(detection_rate_aggr[:,0]==40000)
        aaa_aggr = detection_rate_aggr[ind_val,1][0]
        bbb_aggr = detection_rate_aggr[ind_val,2][0]
        aaa_cons = detection_rate_cons[ind_val,1][0]
        bbb_cons = detection_rate_cons[ind_val,2][0]
        pl.plot(aaa_cons,bbb_cons,label='Conservative (Cone)',linewidth=3,linestyle='--')
        pl.plot(aaa_aggr,bbb_aggr,label='Aggressive (Cone)',linewidth=3,linestyle='-')

    pl.title('D = '+str(int(Dd[ind2]/1e3))+' km')
    pl.xlabel('Slope (deg)')
    #pl.ylabel(r'$\frac{Ntrig(<128antennas)}{Nsimulated showers(135)}$',)
    #pl.ylabel('Ntrig(<128antennas) / Nsimu showers(=135)',)
    pl.ylabel('Trigger rate')
    plt.legend(loc='best',fontsize=10,ncol=2)
    #pl.legend(bbox_to_anchor=(0., -0.3, 0.75, .102), loc=2, ncol=2, mode="expand", borderaxespad=0.,fontsize=10)

if DISPLAY_2D==1:
    dar=(alfa2[-1]-alfa2[0])/(Dd2[-1]-Dd2[0])  *1 #7E-4
    ax = plt.subplot(223)
    #fig1.subplots_adjust(bottom=0.07, hspace=0.3)
    pl.title('Trigger rate - Conservative case')
    pl.ylabel('Distance [km]')
    pl.xlabel('Slope [deg]')

    im = NonUniformImage(ax,interpolation='nearest', extent=(Dd2[0],Dd2[-1],alfa2[0],alfa2[-1]),cmap='jet')
    im.set_data(alfa,Dd,Nconfig_ew_cons)
    ax.images.append(im)
    ax.set_ylim(Dd2[0],Dd2[-1])
    ax.set_xlim(alfa2[0],alfa2[-1])
    im.set_clim(0,np.max(Nconfig_ew_aggr))
    pl.colorbar(im)
    #pl.colorbar(im,orientation="horizontal")
    #pl.gca().set_aspect(dar,adjustable='box')

    ax = plt.subplot(224)
    #fig1.subplots_adjust(bottom=0.07, hspace=0.3)
    pl.title('Trigger rate - Aggressive case')
    pl.ylabel('Distance [km]')
    pl.xlabel('Slope [deg]')

    im = NonUniformImage(ax,interpolation='nearest', extent=(Dd2[0],Dd2[-1],alfa2[0],alfa2[-1]),cmap='jet')
    im.set_data(alfa,Dd,Nconfig_ew_aggr)
    ax.images.append(im)
    ax.set_ylim(Dd2[0],Dd2[-1])
    ax.set_xlim(alfa2[0],alfa2[-1])
    im.set_clim(0,np.max(Nconfig_ew_aggr))
    pl.colorbar(im)
    #pl.colorbar(im,orientation="horizontal")
    #pl.gca().set_aspect(dar,adjustable='box')

    #The following are for different type of plotting methods with other interpolation or smoothing
    '''
    im = pl.imshow((np.transpose(Nconfig_ew_aggr)),cmap='jet',extent=[Dd2[0],Dd2[-1],alfa2[0],alfa2[-1]],origin='lowerleft',aspect=dar)
    im.set_interpolation('spline36')
    #im.set_interpolation('nearest')
    pl.colorbar(orientation="horizontal")
    '''

    '''
    ax = pl.subplot(313) #,adjustable='box',aspect=1)
    pl.title('Trigger rate')
    pl.xlabel('Distance [km]')
    pl.ylabel('Slope [deg]')
    Dd_grid,alfa_grid = np.meshgrid(Dd,alfa)

    #pp = pl.pcolor(Dd_grid,alfa_grid,np.transpose(Nconfig_ew_aggr),cmap='jet')
    pp = pl.pcolormesh(Dd_grid,alfa_grid,np.transpose(Nconfig_tot_cons),cmap='jet',shading='gouraud')
    pl.xlim(Dd2[0],Dd2[-1])
    pl.ylim(alfa2[0],alfa2[-1])
    pl.colorbar(orientation="horizontal")
    pl.gca().set_aspect(dar,adjustable='box')
    '''

    figname = save_folder+'/'+'trigger-rate_'+EE+'eV'+suffix+'.png'
    pl.savefig(figname,dpi=750)
    pl.show()

print 'Finished'

####################################################################################
####################################################################################
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
