#!/usr/bin/env python

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
import danton
import StringIO
import random
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.stats as stats

####################################################################################
# Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg]).
GEOMAGNET = (56.5, 63.18, 2.72)
AZIMUTH = 180.0 #North direction origin
CORE=(0.0,0.0)
TASK_INDEX = 0
part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
nbin=50 #number of bins in the histograms

##########################################################################################################
def main():
    """ Main script allowing to produce the ZHAireS input files for each one of the showers from a DANTON library """

    ##################################################################
    # Test arguments
    if (len(sys.argv)<4 or len(sys.argv)>4):
        print """\
    This script will read the danton put txt file, perform a random selection in the full shower list to get
    a subset of it. Then the script will plot the injection height, elevation and shower energy distribution
    with both the full and selected population along with the fit of the distribution with a gamma function.

    Inputs :
        EE = Neutrino energy e in GeV
        Ntry = number of random draw one want to do among them one will select the subsample of showers
        Nsh = number of shower one want to select

    Ouput:
        The script will produce a plot (3 subplots) and a txt file containing the sublist of selected shower showerID.

    Usage:  python read_danton_output.py EE Ntry Nsh

    Notes:
        The global variable DISPLAY allows you to turn on/off the display of the plots.
    """
        sys.exit(1)

    ####################################################################################
    #Input parameters
    EE=str(sys.argv[1]) #Neutrino energy e
    Ntry = int(sys.argv[2])+1 #number of random draw one want to do among them one will select the subsample of showers
    Nsh = int(sys.argv[3])+1 #number of shower one want to select
    DISPLAY=0

    #Output dir
    ROOT = "/Users/nrenault/Desktop/GRAND/ToyModel/"
    DBDIR =ROOT+"Danton_database/"   #output of danton
    FIGDIR = ROOT+"shower_selection/""
    file=DBDIR+'decays.'+EE+'_all_filtered.txt' #'_all.txt' #'.txt'

    # Parse the events and build up some statistics.
    data=[]
    r1=[]
    r2=[]
    R2=[]
    u2=[]
    for event in danton.iter_event(file):
            lastid = event.id
            for decay in event.decay :
                    r1 = decay.tau_i.position
                    depth = danton.EARTH_RADIUS - np.linalg.norm(r1)
                    r2 = decay.tau_f.position
                    R2 = np.linalg.norm(r2)
                    height = R2 - danton.EARTH_RADIUS
                    u2 = decay.tau_f.direction
                    theta = np.degrees(np.pi - np.arccos(np.dot(u2, r2) / R2)) #in Zhaires convention
                    c = np.dot(u2, event.primary.direction)
                    if c > 1: delta = 0.
                    else: delta = np.arccos(c)
                    dataprod = []
                    et=0.0
                    for i in decay.product:
                        idp = i[0]
                        pp=i[1]
                        ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2) #in GeV ignoring mass energy
                        up=pp/ep
                        thetap=np.degrees(np.pi - np.arccos(np.dot(up, r2) / R2))
                        dataprod.append((event.id,idp,up[0],up[1],up[2],ep,thetap,theta,height,depth,event.primary.energy,decay.tau_f.energy))
                        et=et+ep
                    et=et*1e9 #eV
                    dataprod = np.array(dataprod)

                    data.append((event.id, decay.tau_f.energy, event.primary.energy, depth,
                                 height, theta, delta, decay.generation,et, et/decay.tau_f.energy))

    data = np.array(data)
    print "+ {:} tau decays for {:} incoming neutrinos".format(len(data), lastid+1)
    print('******************************')

    #compute the fraction of earth-skimming/atmospheric neutrino showers
    def compute_prob(K0, K1,dat):
            """Compute the probablity for an event, e.g. conversion in ground.
            """
            r0 = sum(dat[K0,-1])
            r1 = sum(dat[K1,-1])
            s = 1./(r0 + r1)
            return r0 * s, r1 * s

    Kg = data[:,3] > 0.  # Conversion in ground
    Ka = data[:,3] <= 0. # Conversion in the atmosphere
    rg, ra = compute_prob(Kg, Ka,data)
    print "+ {:.0f} % of the conversions are in air".format(100. * ra)

    # Select the earth-skimming neutrinos
    ind = data[:,3]>0
    data_skim = np.array(data[ind,:])

    # Create the data distribution histograms
    hist_theta, bin_edges_theta = np.histogram(data_skim[:,5],bins=nbin,density=True)
    bin_centres_theta = (bin_edges_theta[:-1] + bin_edges_theta[1:])/2
    hist_injh, bin_edges_injh = np.histogram(np.log10(data_skim[:,4]),bins=nbin,density=True)
    bin_centres_injh = (bin_edges_injh[:-1] + bin_edges_injh[1:])/2
    hist_et, bin_edges_et = np.histogram(np.log10(data_skim[:,8]),bins=nbin,density=True)
    bin_centres_et = (bin_edges_et[:-1] + bin_edges_et[1:])/2

    # Define model function to be used to fit to the data above:
    def poisson(k, *p):
        lamb,shift = p
        return (lamb**(k-shift)/factorial((k-shift))) * np.exp(-lamb)
    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

    # Fit those distributions
    #coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=[0.18,90., 1.])
    #coeff, var_matrix = curve_fit(poisson, bin_centres, hist,p0=[1.,90.])
    fit_alpha_theta, fit_loc_theta, fit_beta_theta=stats.gamma.fit(data_skim[:,5])
    fit_alpha_injh, fit_loc_injh, fit_beta_injh=stats.gamma.fit(np.log10(data_skim[:,4]))
    fit_alpha_et, fit_loc_et, fit_beta_et=stats.gamma.fit(np.log10(data_skim[:,8]))

    # Get the fitted curve
    hist_fit_theta = stats.gamma.pdf(bin_centres_theta,fit_alpha_theta, loc=fit_loc_theta, scale=fit_beta_theta)
    hist_fit_injh = stats.gamma.pdf(bin_centres_injh,fit_alpha_injh, loc=fit_loc_injh, scale=fit_beta_injh)
    hist_fit_et = stats.gamma.pdf(bin_centres_et,fit_alpha_et, loc=fit_loc_et, scale=fit_beta_et)

    #Filter showers among which the selection will be performed
    ind2 = np.logical_and(data_skim[:,5]-90.<=15.,data_skim[:,5]-90.>=0.) #np.array(ind)
    data_skim2 = np.array(data_skim[ind2,:])
    ind_run1 = np.logical_and(data_skim2[:,4]>=0.,data_skim2[:,4]<=13000.) #np.array(ind)
    data_inj = np.array(data_skim2[ind_run1,:])

    #Set a list of showerID, you want to see included in the final list of selected showers
    list_run2 = np.array([],dtype=float) #np.array([1000689,1002953,1003115],dtype=float)
    ix2 = np.in1d(data_inj[:,0], list_run2)
    id_run2 = np.where(ix2)[0]
    id_run = np.copy(id_run2)

    #Perform the random draw Ntry times and register the shower list everytime some criteria are met
    theta_to_select = np.copy(data_inj[:,5])
    theta_to_select = np.delete(theta_to_select,id_run2)
    proba_theta = stats.gamma.pdf(theta_to_select,fit_alpha_theta, loc=fit_loc_theta, scale=fit_beta_theta)
    for iit in range(0,Ntry):
        #Perform the draw
        id_rdm_sel = np.random.choice(range(0, len(theta_to_select)),replace=False,size=Nsh-len(id_run)) #,p=proba_theta/np.sum(proba_theta)) #np.ones(np.shape(proba_theta))/np.shape(proba_theta)) #proba_theta/np.sum(proba_theta))
        id_rdm_sel = np.concatenate((id_run,id_rdm_sel),axis=0)
        data_studied = np.array(data_inj[id_rdm_sel,:])
        results = map(int, np.reshape(data_studied[:,0],(np.shape(data_studied)[0])))

        #Test if this selection meets one/some criterion/a
        if np.sum(np.array(pg_theta[0][1])>0.30)==0: #0.35
            #Eventually plots the corresponding distribution
            fig1 = pl.figure(figsize=(5*3.13,3.9*3.13))
            pl.subplot(221)
            pg_theta = pl.hist([data_skim[:,5],data_studied[:,5]],bins=nbin,normed=1,label=['Full (histo)','Selected'])
            pl.plot(bin_centres_theta, hist_theta,'-k', label='Full (line)',linewidth=2)
            pl.plot(bin_centres_theta, hist_fit_theta,'--r', label='Fit',linewidth=2)
            pl.xlabel('Elevation')
            pl.ylabel('#')
            pl.grid(True)
            pl.legend()

            pl.subplot(222)
            pg_injh = pl.hist([np.log10(data_skim[:,4]),np.log10(data_studied[:,4])],bins=nbin,normed=1,label=['Full (histo)','Selected'])
            pl.plot(bin_centres_injh, hist_injh,'-k', label='Full (line)',linewidth=2)
            pl.plot(bin_centres_injh, hist_fit_injh,'--r', label='Fit',linewidth=2)
            pl.xlabel('Injection height')
            pl.ylabel('#')
            pl.grid(True)
            pl.legend()

            pl.subplot(223)
            pg_et = pl.hist([np.log10(data_skim[:,8]),np.log10(data_studied[:,8])],bins=nbin,normed=1,label=['Full (histo)','Selected'])
            pl.plot(bin_centres_et, hist_et,'-k', label='Full (line)',linewidth=2)
            pl.plot(bin_centres_et, hist_fit_et,'--r', label='Fit',linewidth=2)
            pl.xlabel('Shower energy')
            pl.ylabel('#')
            pl.grid(True)
            pl.legend()

            #Save the distributions and shower list
            if not os.path.exists(FIGDIR):
                os.makedirs(FIGDIR)
            figname = FIGDIR+'plots_EE'+EE+'_sel'+str(iit)+'.png'
            pl.savefig(figname,dpi=350)
            showers_list_file = ROOT+'/shower_selection/liste_EE'+EE+'_sel'+str(iit)+'.txt'
            np.savetxt(showers_list_file,results,fmt='%d' ,header="ShowerID")
            if DISPLAY==0:
                pl.close(fig1)
            elif DISPLAY==1:
                pl.show()

##########################################################################################################
if __name__ == '__main__':
    main()
