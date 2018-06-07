# Analyse ZHaireS + ComputeVoltage outputs
# Produced by Nicolas (see email on 24/10/2017)

import sys
import numpy as np
import pylab as pl
pl.ioff()

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<2 or len(sys.argv)>2):
    print """\
        This script will either plots the distributions of the generated and simulated CR populations, or compute and plot the performances of GRANDproto300 for the CR detection.

        Usage:  python CR_calc_GP300.py [opt: distrib/ for CR distributions or perf for GRANDproto300 performances]
        """
    sys.exit(1)
###########################

###########################################################################
# Set some parameters
sec2yr = 1./(3600.*24.*365.25)
pi = 3.1415927
nantmin = 5 #8
ndays = 1. #in days
dt_day = 3600.*24.*ndays #in sec
#elow = np.array([pow(10,float(x)/10) for x in range(160,int(np.ceil(np.log10(1e17)*10)))])/1e9
#ehigh = np.array([pow(10,float(x)/10) for x in range(int(np.ceil(np.log10(1e17)*10)),200)])/1e9

J1 = 1.06e29  # /s/m2/sr/eV
gamma1 = 3.26
J2 = 4.5e17
gamma2 = 2.65
eth=np.array([pow(10,float(x)/10) for x in range(160,200)])
eankle = [pow(10,18.26)]
leth = np.concatenate((eth[np.where(eth<=eankle)],eankle))
heth = np.concatenate((eankle,eth[np.where(eth>=eankle)]))
f1 = [J1*pow(x,-gamma1)*pow(x,3)/1e24 for x in leth]  # below ankle
f2 = [J2*pow(x,-gamma2)*pow(x,3)/1e24 for x in heth] # above ankle

###########################################################################
# Load results from txt files
#Set1
datadir = '/Users/nrenault/Desktop/GRAND/CRs_GP300/simus_set1/detection/'
resfile = datadir+'detection_count_50-200MHz.txt'
res1 = np.loadtxt(resfile)
drawfile = '/Users/nrenault/Desktop/GRAND/CRs_GP300/inp_set1/summary_table.txt'
draw1 = np.loadtxt(drawfile)
#Set2
datadir = '/Users/nrenault/Desktop/GRAND/CRs_GP300/simus_set2/detection/'
resfile = datadir+'detection_count_50-200MHz.txt'
res2 = np.loadtxt(resfile)
drawfile = '/Users/nrenault/Desktop/GRAND/CRs_GP300/inp_set2/summary_table.txt'
draw2 = np.loadtxt(drawfile)
#Set4
datadir = '/Users/nrenault/Desktop/GRAND/CRs_GP300/simus_set3/detection/'
resfile = datadir+'detection_count_50-200MHz.txt'
res3 = np.loadtxt(resfile)
drawfile = '/Users/nrenault/Desktop/GRAND/CRs_GP300/inp_set3/summary_table.txt'
draw3 = np.loadtxt(drawfile)
#Set4
datadir = '/Users/nrenault/Desktop/GRAND/CRs_GP300/simus_set4/detection/'
resfile = datadir+'detection_count_50-200MHz.txt'
res4 = np.loadtxt(resfile)
drawfile = '/Users/nrenault/Desktop/GRAND/CRs_GP300/inp_set4/summary_table.txt'
draw4 = np.loadtxt(drawfile)

###########################################################################
#Correct the evtID to account for the multiple txt files
dimTry1 = int(np.amax(draw1[:,1]))
draw2[np.where(draw2[:,1]!=0)[0],1] = draw2[np.where(draw2[:,1]!=0)[0],1]+dimTry1
res2[:,1] = res2[:,1]+dimTry1
dimTry2 = int(np.amax(draw2[:,1]))
draw3[np.where(draw3[:,1]!=0)[0],1] = draw3[np.where(draw3[:,1]!=0)[0],1]+dimTry2
res3[:,1] = res3[:,1]+dimTry2
#We do not do the same for the set4 because set4 covers a different energy range
draw_tmp = np.concatenate((draw1,draw2,draw3,draw4))
res_tmp = np.concatenate((res1,res2,res3,res4))

###########################################################################
#Sort the array with increasing energy while keeping the order of the events
eny = np.sort(np.unique(draw_tmp[:,2]))
draw = np.zeros(np.shape(draw_tmp))
istart = 0
for ie in range(len(eny)):
    ind = np.where(draw_tmp[:,2]==eny[ie])[0]#[0]
    draw[istart:istart+np.size(ind),:] = draw_tmp[ind,:]
    istart = istart+np.size(ind)

eny = np.sort(np.unique(res_tmp[:,0]))
res = np.zeros(np.shape(res_tmp))
istart = 0
for ie in range(len(eny)):
    ind = np.where(res_tmp[:,0]==eny[ie])[0]#[0]
    res[istart:istart+np.size(ind),:] = res_tmp[ind,:]
    istart = istart+np.size(ind)

###########################################################################
#Load infos from draw file
Nevt = draw[:,0]
num_selevt = draw[:,1] #Ntry = np.array(draw[:,4],dtype=int)
logEny = np.array(np.log10(draw[:,2]*1e18),dtype=float) #EeV
azim = draw[:,3] #GRAND conv
zen = draw[:,4] #GRAND conv
Nant_footprint = draw[:,5]
#xmean = draw[:,6]
#ymean = draw[:,7]
xc = draw[:,8]
yc = draw[:,9]
altc = draw[:,10]
dx = draw[:,11]
dy = draw[:,12]
#GdAlt = draw[:,13]
#MinAnt = draw[:,14]

ind_sel = np.where(num_selevt!=0)[0] #selected for simus, From summary tab file
Ntry = np.zeros(len(ind_sel),dtype=int)
Ntry[0] = Nevt[ind_sel[0]]
Ntry[1:] = [Nevt[ind_sel[i]]-Nevt[ind_sel[i-1]] for i in range(1,len(ind_sel))]
dimTry = int(np.amax(num_selevt))
dimE = len(ind_sel)/dimTry
Ntry = Ntry.reshape((dimE,dimTry)).T

#Nevt_sel = np.array(Nevt[ind_sel],dtype=int) #
ID_selevt = np.array(num_selevt[ind_sel],dtype=int).reshape((dimE,dimTry)).T #ID selected for simus, From summary tab file
phi = np.array(azim[ind_sel],dtype=float).reshape((dimE,dimTry)).T
theta = np.array(180.-zen[ind_sel],dtype=float).reshape((dimE,dimTry)).T
Elog = np.around(logEny[ind_sel],1).reshape((dimE,dimTry)).T
eny = np.sort(np.unique(Elog))
ee = pow(10,np.unique(Elog))
ind_le = np.where(ee<=eankle)
ind_he = np.where(ee>eankle)
Adraw = 2*dx[0]/1e3*2*dy[0]/1e3  # Try area (km2)
surf = 175 #km2 #175 antennas with a 1km-step, 250m- and 500m-step arrays are included within it
Aproj = np.cos(theta*pi/180)*Adraw

###########################################################################
# Load info for detection_count file
evt_number_sim = np.array(res[:,1],dtype=int) #ID selected for simus, From detection count file
eny_sim = res[:,0] #From detection count file #to be compared to E
Nant_ew_cons = res[:,2]
Nant_ns_cons = res[:,3]
Nant_up_cons = res[:,4]
Nant_tot_cons = res[:,5]
Nant_ew_agg = res[:,6]
Nant_ns_agg = res[:,7]
Nant_up_agg = res[:,8]
Nant_tot_agg = res[:,9]

########################################################################################################
########################################################################################################
########################################################################################################
def run():
    #print np.unique(evt_number_sim)
    #print np.unique(ID_selevt)
    #stop
    #Compute CR Aeff,aperture,exposure
    aper_cons = np.zeros(dimE)
    aper_agg = np.zeros(dimE)
    aper_opt = np.zeros(dimE)
    expo_cons = np.zeros(dimE)
    expo_agg = np.zeros(dimE)
    expo_opt = np.zeros(dimE)

    Aeff_phi_cons = np.zeros((dimE,dimTry))
    Aeff_phi_agg = np.zeros((dimE,dimTry))
    Aeff_phi_opt = np.zeros((dimE,dimTry))
    Aeff_theta_phi_opt = np.zeros((dimE,6)) #6=dim_th related to the size of the hard coded th_bin in fction compute_Aeff_theta
    Aeff_theta_phi_agg = np.zeros((dimE,6)) #6=dim_th related to the size of the hard coded th_bin in fction compute_Aeff_theta
    Aeff_theta_phi_cons = np.zeros((dimE,6)) #6=dim_th related to the size of the hard coded th_bin in fction compute_Aeff_theta
    for ie in range(dimE):
        th = np.sort(theta[:,ie])
        print 'E=',Elog[0,ie]
        for ievt in range(dimTry): #loop over theta usually but here each event has its own theta so it is the same
            sel = np.where(np.logical_and(eny_sim == Elog[ievt,ie], evt_number_sim==ievt+1))[0] #sel = np.intersect1d(np.where(eny_sim == eu[ie]),np.where(evt_number_sim == ievt))
            ok_agg = np.logical_or(Nant_ew_agg[sel]>=nantmin,Nant_ns_agg[sel]>=nantmin)#,Nant_up_agg[0][sel]>=nantmin)
            ok_cons = np.logical_or(Nant_ew_cons[sel]>=nantmin,Nant_ns_cons[sel]>=nantmin)#,Nant_up_cons[0][sel]>=nantmin)
            Ragg = float(np.sum(ok_agg))/np.sum(Ntry[ievt,ie])
            Rcons = float(np.sum(ok_cons))/np.sum(Ntry[ievt,ie])

            ith = np.where(th==theta[ievt,ie])[0][0]
            ieny = np.where(eny==Elog[ievt,ie])[0][0]
            Aeff_phi_opt[ieny,ith] = Aproj[ievt,ie]
            Aeff_phi_agg[ieny,ith] = Aeff_phi_opt[ieny,ith]*Ragg   # Effective area integrated over phi
            Aeff_phi_cons[ieny,ith] = Aeff_phi_opt[ieny,ith]*Rcons   # Effective area integrated over phi

        th_bin,Aeff_theta_phi_opt[ie,:],Aeff_theta_phi_agg[ie,:],Aeff_theta_phi_cons[ie,:] = compute_Aeff_theta(th,Aeff_phi_opt[ie,:],Aeff_phi_agg[ie,:],Aeff_phi_cons[ie,:])

        #Idea : theta must be sorted from low to high value and not in disorder as it is now so use th or similar
        aper_agg[ie] = 2*pi*np.trapz(Aeff_phi_agg[ie,:]*np.sin(th*pi/180),th*pi/180)
        aper_cons[ie] = 2*pi*np.trapz(Aeff_phi_cons[ie,:]*np.sin(th*pi/180),th*pi/180)
        aper_opt[ie] = 2*pi*np.trapz(Aeff_phi_opt[ie,:]*np.sin(th*pi/180),th*pi/180)

        expo_agg[ie] = aper_agg[ie]*dt_day
        expo_cons[ie] = aper_cons[ie]*dt_day
        expo_opt[ie] = aper_opt[ie]*dt_day

    # 1e6 to transfer from km2 to m2
    evt1d_agg = 1e6*(np.trapz(expo_agg[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_agg[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
    evt1d_cons = 1e6*(np.trapz(expo_cons[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_cons[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
    evt1d_opt = 1e6*(np.trapz(expo_opt[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_opt[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
    evt1d_agg_diff = 1e9*1e6*np.concatenate([expo_agg[ind_le]*J1*pow(ee[ind_le],-gamma1),expo_agg[ind_he]*J2*pow(ee[ind_he],-gamma2)])  #GeV-1 day-1
    evt1d_cons_diff = 1e9*1e6*np.concatenate([expo_cons[ind_le]*J1*pow(ee[ind_le],-gamma1),expo_cons[ind_he]*J2*pow(ee[ind_he],-gamma2)])

    evt1d_agg_1e19 = 1e6*np.trapz(expo_agg[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])
    evt1d_cons_1e19 = 1e6*np.trapz(expo_cons[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])
    evt1d_opt_1e19 = 1e6*np.trapz(expo_opt[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])

    evt1d_agg_1e1819 = 1e6*np.trapz(expo_agg[3:5]*J2*pow(ee[3:5],-gamma2),ee[3:5])
    evt1d_cons_1e1819 = 1e6*np.trapz(expo_cons[3:5]*J2*pow(ee[3:5],-gamma2),ee[3:5])
    evt1d_opt_1e1819 = 1e6*np.trapz(expo_opt[3:5]*J2*pow(ee[3:5],-gamma2),ee[3:5])

    evt1d_agg_lastbd = 1e6*np.trapz(np.array([expo_agg[-1],expo_agg[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))
    evt1d_cons_lastbd = 1e6*np.trapz(np.array([expo_cons[-1],expo_cons[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))
    evt1d_opt_lastbd = 1e6*np.trapz(np.array([expo_opt[-1],expo_opt[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))

    evt1d_agg_1e16 = 1e6*np.trapz(np.array([expo_agg[0],expo_agg[0]])*J1*pow(np.array([10**16.25,10**16.75]),-gamma1),np.array([10**16.25,10**16.75]))
    evt1d_cons_1e16 = 1e6*np.trapz(np.array([expo_cons[0],expo_cons[0]])*J1*pow(np.array([10**16.25,10**16.75]),-gamma1),np.array([10**16.25,10**16.75]))
    evt1d_opt_1e16 = 1e6*np.trapz(np.array([expo_opt[0],expo_opt[0]])*J1*pow(np.array([10**16.25,10**16.75]),-gamma1),np.array([10**16.25,10**16.75]))

    evt1d_agg_1e17 = 1e6*np.trapz(np.array([expo_agg[1],expo_agg[1]])*J1*pow(np.array([10**16.75,10**17.25]),-gamma1),np.array([10**16.75,10**17.25]))
    evt1d_cons_1e17 = 1e6*np.trapz(np.array([expo_cons[1],expo_cons[1]])*J1*pow(np.array([10**16.75,10**17.25]),-gamma1),np.array([10**16.75,10**17.25]))
    evt1d_opt_1e17 = 1e6*np.trapz(np.array([expo_opt[1],expo_opt[1]])*J1*pow(np.array([10**16.75,10**17.25]),-gamma1),np.array([10**16.75,10**17.25]))

    evt1d_agg_1e175 = 1e6*np.trapz(np.array([expo_agg[2],expo_agg[2]])*J1*pow(np.array([10**17.25,10**17.75]),-gamma1),np.array([10**17.25,10**17.75]))
    evt1d_cons_1e175 = 1e6*np.trapz(np.array([expo_cons[2],expo_cons[2]])*J1*pow(np.array([10**17.25,10**17.75]),-gamma1),np.array([10**17.25,10**17.75]))
    evt1d_opt_1e175 = 1e6*np.trapz(np.array([expo_opt[2],expo_opt[2]])*J1*pow(np.array([10**17.25,10**17.75]),-gamma1),np.array([10**17.25,10**17.75]))

    evt1d_agg_1e18 = 1e6*np.trapz(np.array([expo_agg[3],expo_agg[3]])*J1*pow(np.array([10**17.75,eankle[0]]),-gamma1),np.array([10**17.75,eankle[0]]))
    evt1d_cons_1e18 = 1e6*np.trapz(np.array([expo_cons[3],expo_cons[3]])*J1*pow(np.array([10**17.75,eankle[0]]),-gamma1),np.array([10**17.75,eankle[0]]))
    evt1d_opt_1e18 = 1e6*np.trapz(np.array([expo_opt[3],expo_opt[3]])*J1*pow(np.array([10**17.75,eankle[0]]),-gamma1),np.array([10**17.75,eankle[0]]))

    evt1d_agg_eankle = 1e6*np.trapz(np.array([expo_agg[4],expo_agg[4]])*J2*pow(np.array([eankle[0],10**18.75]),-gamma2),np.array([eankle[0],10**18.75]))+1e6*np.trapz(np.array([expo_agg[4],expo_agg[4]])*J2*pow(np.array([10**18.75,10**19.25]),-gamma2),np.array([10**18.75,10**19.25]))+1e6*np.trapz(np.array([expo_agg[5],expo_agg[5]])*J2*pow(np.array([10**19.25,10**19.75]),-gamma2),np.array([10**19.25,10**19.75]))
    evt1d_cons_eankle = 1e6*np.trapz(np.array([expo_cons[4],expo_cons[4]])*J2*pow(np.array([eankle[0],10**18.75]),-gamma2),np.array([eankle[0],10**18.75]))+1e6*np.trapz(np.array([expo_cons[4],expo_cons[4]])*J2*pow(np.array([10**18.75,10**19.25]),-gamma2),np.array([10**18.75,10**19.25]))+1e6*np.trapz(np.array([expo_cons[5],expo_cons[5]])*J2*pow(np.array([10**19.25,10**19.75]),-gamma2),np.array([10**19.25,10**19.75]))
    evt1d_opt_eankle = 1e6*np.trapz(np.array([expo_opt[4],expo_opt[4]])*J2*pow(np.array([eankle[0],10**18.75]),-gamma2),np.array([eankle[0],10**18.75]))+1e6*np.trapz(np.array([expo_opt[4],expo_opt[4]])*J2*pow(np.array([10**18.75,10**19.25]),-gamma2),np.array([10**18.75,10**19.25]))+1e6*np.trapz(np.array([expo_opt[5],expo_opt[5]])*J2*pow(np.array([10**19.25,10**19.75]),-gamma2),np.array([10**19.25,10**19.75]))

    ndays_ref = 365.*1.  # duration of observation (days)
    dt_ref = 3600.*24.*ndays_ref #in s
    rescale_dt = dt_ref/dt_day
    expo_agg = expo_agg *rescale_dt# *sec2yr
    expo_cons = expo_cons *rescale_dt# *sec2yr
    expo_opt = expo_opt *rescale_dt# *sec2yr

    print 'expo_cons=',expo_cons
    print 'expo_agg=',expo_agg
    print 'expo_opt=',expo_opt

    #print 'Aperture (agressive):',aper_agg,'km2.sr   -   Exposure (aggressive):',expo_agg*sec2yr,'km2.sr.yr   -   Event rate (agressive):',evt_eny_agg
    #print 'Aperture (conservative):',aper_cons,'km2.sr   -   Exposure (conservative):',expo_cons*sec2yr,'km2.sr.yr   -   Event rate (conservative):',evt_eny_cons
    #print 'Optimal apperture:',aper_opt,'km2.sr   -   Exposure (optimal):',expo_opt*sec2yr,'km2.sr.yr   -   Event rate (optimal):',evt_eny_opt
    #print ' '

    #print 'exposure_300 = ',expo_cons[-1]*np.arange(1,12) #*365.*3600.*24./dt_day  *sec2yr

    print 'Expected yearly event rate in 10^17-10^19.5eV (agressive):',evt1d_agg
    print 'Expected yearly event rate in 10^17-10^19.5eV (conservative):',evt1d_cons
    print ' '
    print 'Expected yearly event rate in [10^18,10^19[ eV (agressive):',evt1d_agg_1e1819
    print 'Expected yearly event rate in [10^18,10^19[ eV (conservative):',evt1d_cons_1e1819
    print ' '
    print 'Expected yearly event rate for E>=10^19eV (agressive):',evt1d_agg_1e19
    print 'Expected yearly event rate for E>=10^19eV (conservative):',evt1d_cons_1e19
    print ' '
    print 'Expected yearly event rate for E>=10^19.5eV (agressive):',evt1d_agg_lastbd
    print 'Expected yearly event rate for E>=10^19.5eV (conservative):',evt1d_cons_lastbd
    print ' '
    print 'Expected yearly event rate around 10^16 eV (agressive):',evt1d_agg_1e16 *rescale_dt
    print 'Expected yearly event rate around 10^16 eV (conservative):',evt1d_cons_1e16 *rescale_dt
    print ' '
    print 'Expected yearly event rate around 10^17 eV (agressive):',evt1d_agg_1e17 *rescale_dt
    print 'Expected yearly event rate around 10^17 eV (conservative):',evt1d_cons_1e17 *rescale_dt
    print ' '
    print 'Expected yearly event rate around 10^17.5 eV (agressive):',evt1d_agg_1e175 *rescale_dt
    print 'Expected yearly event rate around 10^17.5 eV (conservative):',evt1d_cons_1e175 *rescale_dt
    print ' '
    print 'Expected yearly event rate around 10^18eV (agressive):',evt1d_agg_1e18 *rescale_dt
    print 'Expected yearly event rate around 10^18eV (conservative):',evt1d_cons_1e18 *rescale_dt
    print ' '
    print 'Expected yearly event rate for E>=10^18.26 eV (agressive):',evt1d_agg_eankle *rescale_dt
    print 'Expected yearly event rate for E>=10^18.26 eV (conservative):',evt1d_cons_eankle *rescale_dt

    expo_auger = np.array([5400*dt_day*rescale_dt]) #np.array([48000./(9.*365.*3600.*24.),48000./(9.*365.*3600.*24.)])
    evt1d_cons_auger = 1e6*np.trapz(expo_auger*J2*pow(np.array([ee[-1:],1e20]),-gamma2),np.array([ee[-1:],1e20]))

    ###############################################################
    ##### Plots

    # Observed Flux
    pl.figure(17)
    pl.loglog(leth,f1,lw=3,label='flux below ankle')
    pl.loglog(heth,f2,lw=3,label='flux above ankle')
    pl.xlabel('Energy (eV)')
    pl.ylabel('Flux*E$^3$ (eV$^2$/s/sr/m$^2$)')
    pl.legend(loc='best')
    pl.title("TA flux - astro-ph 1511.07510")
    pl.grid(True)

    # Aperture, exposure
    pl.figure(1)
    pl.loglog(ee,aper_agg,'+-',lw=3,label='Agressive',mew=3, ms=10)
    pl.loglog(ee,aper_cons,'+-',lw=3,label='Conservative',mew=3, ms=10)
    pl.loglog([1e19,1e20],[5400,5400],'-',lw=3,label='Auger',mew=3,ms=10)
    pl.xlabel('Energy (eV)')
    pl.ylabel('Aperture (km$^2$.sr)')
    pl.grid(True)
    pl.legend(loc='best')
    pl.title("GRANDproto300 CR aperture")

    pl.figure(2)
    pl.loglog(ee,expo_agg,'+-',lw=3,label='Agressive',mew=3, ms=10)
    pl.loglog(ee,expo_cons,'+-',lw=3,label='Conservative',mew=3, ms=10)
    pl.loglog([1e19,1e20],[48000/sec2yr,48000/sec2yr],'-',lw=3,label='Auger 9yr',mew=3,ms=10)
    pl.xlabel('Energy (eV)')
    pl.ylabel('Exposure (km$^2$.sr.s)')
    pl.grid(True)
    pl.legend(loc='best')
    pl.title("GRANDProto300 CR exposure")

    # Effective area
    pl.subplots(nrows=1, ncols=2,figsize=(5.*3.13,3.1*3.13))
    ax = pl.subplot(121)
    thc = (th_bin[1:]+th_bin[0:-1])/2.
    pl.semilogy(thc,np.cos(thc*pi/180)*Adraw,'-',lw=3)#,label='Optimal')
    #pl.semilogy(thc,Aeff_theta_phi_opt[3,:],'x-',lw=3,label='Optimal')
    for ie in range(dimE):
        pl.semilogy(thc,Aeff_theta_phi_agg[ie,:],'x-',lw=3,mew=3, ms=10)#,label='Aggressive - '+str(np.round(10**Elog[0,ie]/1e18,1))+' EeV')
    pl.xlabel('Theta (deg)')
    pl.ylabel('$<A_{eff}>_\phi$ (km$^2$)')
    pl.grid(True)
    pl.title("GRANDproto300 CR $A_{eff}$ (aggressive case)")
    #pl.legend(loc='lower left',fontsize=11)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])

    ax = pl.subplot(122)
    pl.semilogy(thc,np.cos(thc*pi/180)*Adraw,'-',lw=3,label='Optimal')
    #pl.semilogy(thc,Aeff_theta_phi_opt[3,:],'x-',lw=3,label='Optimal')
    for ie in range(dimE):
        pl.semilogy(thc,Aeff_theta_phi_cons[ie,:],'x-',lw=3,mew=3, ms=10,label=str(np.round(10**Elog[0,ie]/1e18,1))+' EeV') #,label='Conservative - '+str(np.round(10**Elog[0,ie]/1e18,1))+' EeV')
    pl.xlabel('Theta (deg)')
    pl.ylabel('$<A_{eff}>_\phi$ (km$^2$)')
    pl.grid(True)
    pl.title("GRANDproto300 CR $A_{eff}$ (conservative case)")
    #pl.legend(loc='lower left',fontsize=11)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(fontsize=11,loc='upper center', bbox_to_anchor=(-0.1, -0.1),
          fancybox=False, shadow=False, ncol=8)

    # Daily event rate
    pl.figure(5)
    pl.loglog(ee,evt1d_agg_diff,'+-',lw=3,label='Agressive',mew=3, ms=10)
    pl.loglog(ee,evt1d_cons_diff,'+-',lw=3,label='Conservative',mew=3, ms=10)
    pl.xlabel('Energy (eV)')
    pl.ylabel('dN/dE/dt (GeV$^{-1}$ day$^{-1}$)')
    pl.grid(True)
    pl.legend(loc='best')
    pl.title("GRANDProto300 CR event rate")

    pl.show()

########################################################################################################
def compute_Aeff_theta(th,Aeff_phi_opt,Aeff_phi_agg,Aeff_phi_cons):
    th_bin = np.array([60.,65.,70.,75.,80.,85.,90.],dtype=float)
    Aeff_theta_phi_opt = np.zeros(len(th_bin)-1)
    Aeff_theta_phi_agg = np.zeros(len(th_bin)-1)
    Aeff_theta_phi_cons = np.zeros(len(th_bin)-1)
    for ith in range(len(th_bin)-1):
        ind = np.where(np.logical_and(th>=th_bin[ith],th<th_bin[ith+1]))[0]
        if len(ind)>0:
            Aeff_theta_phi_opt[ith]=np.mean(Aeff_phi_opt[ind])
            Aeff_theta_phi_agg[ith]=np.mean(Aeff_phi_agg[ind])
            Aeff_theta_phi_cons[ith]=np.mean(Aeff_phi_cons[ind])

    return th_bin,Aeff_theta_phi_opt,Aeff_theta_phi_agg,Aeff_theta_phi_cons

########################################################################################################
def distrib_events(ind_sel,logeny,azim,zen,xc,yc,altc,Nant_footprint,Ntry):

    # Histo energy
    fig1 = pl.figure(figsize=(21.0,9.7))
    pl.subplot(131)
    pl.hist([logeny,logeny[ind_sel]],bins=[16.25,16.75,17.25,17.75,18.25,18.75,19.25,19.75],density=1,label=['Generated','Simulated'],align='mid')
    pl.xlabel('Energy')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Cosmic-ray shower energy')

    pl.subplot(132)
    pl.hist([azim,azim[ind_sel]],bins=25,density=1,label=['Generated','Simulated'])
    pl.xlabel('GRAND azimuth')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Cosmic-ray shower azimuth in GRAND convention')

    pl.subplot(133)
    pl.hist([zen,zen[ind_sel]],bins=25,density=1,label=['Generated','Simulated'])
    pl.xlabel('GRAND zenith')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Cosmic-ray shower zenith in GRAND convention')

    fig4 = pl.figure(figsize=(21.0,9.7))
    pl.subplot(131)
    pl.hist([xc,xc[ind_sel]],bins=25,density=1,label=['Generated','Simulated'])
    pl.xlabel('Core X ')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Shower core X position')

    pl.subplot(132)
    pl.hist([yc,yc[ind_sel]],bins=25,density=1,label=['Generated','Simulated'])
    pl.xlabel('Core Y')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Shower core Y position')

    pl.subplot(133)
    pl.hist([altc,altc[ind_sel]],bins=25,density=1,label=['Generated','Simulated'])
    pl.xlabel('Core Z')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Shower core altitude')

    fig5 = pl.figure(figsize=(21.0,9.7))
    pl.hist([(Nant_footprint[ind_sel])],bins=25,density=0,label=['Simulated'])
    pl.xlabel('log10 Number of antennas in footprint')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Number of antennas in footprint')

    #'''
    fig6 = pl.figure(figsize=(21.0,9.7))
    pl.hist(np.log10(Ntry),bins=25,density=0,label=[str('{0:.2e}'.format(10**eny[ie]))+' eV' for ie in range(len(eny))])
    pl.xlabel('log10 Number of tries')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Number of tries')
    #'''
    pl.show()

########################################################################################################
########################################################################################################
########################################################################################################
if __name__ == '__main__':
    try:
        DoP = sys.argv[1]
    except:
        DoP = 'perf'
    if DoP=='distrib':
        distrib_events(ind_sel,logEny,azim,zen,xc,yc,altc,Nant_footprint,Ntry)
    elif DoP=='perf':
        run()
