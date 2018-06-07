# Analyse ZHaireS + ComputeVoltage outputs
# Produced by Nicolas (see email on 24/10/2017)

import sys
import numpy as np
import pylab as pl
pl.ioff()

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

# Load results
#datadir = '/home/martineau/GRAND/GRAND/data/CRs/'
datadir = '/Users/nrenault/Desktop/GRAND/CR_v1/simus/detection_BU/'
resfile = datadir+'detection_count.txt'
res = np.loadtxt(resfile)

#drawfile = datadir+'/summary_table_evt0to9.txt'
drawfile = '/Users/nrenault/Desktop/GRAND/CR_v1/inp/summary_table_evt0to9.txt'
draw = np.loadtxt(drawfile)

#Store infos from draws
ind = np.where(draw[:,2]!=180.)

num_evt = draw[ind,0]
eny = draw[ind,1] #EeV
azim = draw[ind,2]
zen = draw[ind,3]
Ntry = np.array(draw[ind,4],dtype=int)
Nant_footprint = draw[ind,5]
xc = draw[ind,6]
yc = draw[ind,7]
dx = draw[ind,8]
dy = draw[ind,9]
step = draw[ind,10]
Nx = draw[ind,11]
Ny = draw[ind,12]
slope = draw[ind,13]
mountain_height = draw[ind,14]
GdAlt = draw[ind,15]
MinAnt = draw[ind,16]

#Store infos from detection counts
ind2 = np.where(res[:,1]!=180.)

Neny = np.array(res[ind2,0],dtype=int)
EE = np.array([17.5,18.0,18.5,19.0,19.5],dtype=float)
Elog = np.array([EE[ie] for ie in Neny],dtype=float)
eu = np.unique(Elog)
E= pow(10,Elog)
ee = pow(10,eu)
ind_le = np.where(ee<=eankle)
ind_he = np.where(ee>eankle)

phi = res[ind2,1]
theta = np.array(180.-res[ind2,2])
th = np.unique(theta)
Adraw = 2*dx[0][0]/1e3*2*dy[0][0]/1e3  # Try area (km2)
surf = Nx[0][0]*Ny[0][0]*step[0][0]/1e3*step[0][0]/1e3
Aproj = np.cos(theta*pi/180)*Adraw
#print theta,Adraw,Aproj
evt_number = res[ind2,3]

#Aproj_200k = np.cos(theta*pi/180) *(np.sqrt(200000.)+100.)**2.
Aproj_200k = np.cos(theta*pi/180)  *(np.sqrt(10000.)+100.)**2. *20.

Nant_ew_cons = res[ind2,4]
Nant_ns_cons = res[ind2,5]
Nant_up_cons = res[ind2,6]
Nant_tot_cons = res[ind2,7]
Nant_ew_agg = res[ind2,8]
Nant_ns_agg = res[ind2,9]
Nant_up_agg = res[ind2,10]
Nant_tot_agg = res[ind2,11]

#Compute CR Aeff,aperture,exposure
aper_cons = np.zeros(np.size(eu))
aper_agg = np.zeros(np.size(eu))
aper_opt = np.zeros(np.size(eu))
expo_cons = np.zeros(np.size(eu))
expo_agg = np.zeros(np.size(eu))
expo_opt = np.zeros(np.size(eu))

Aeff_phi_cons = np.zeros((np.size(eu),np.size(th)))
Aeff_phi_agg = np.zeros((np.size(eu),np.size(th)))
Aeff_phi_opt = np.zeros((np.size(eu),np.size(th)))
Aeff_phi_opt_200k = np.zeros((np.size(eu),np.size(th)))

for i in range(np.size(eu)):
  e = eu[i]
  for t in range(np.size(th)):
    print "theta=",th[t],"E=",e
    #sel = np.where((Elog == e) & (theta == th[t]))
    sel = np.intersect1d(np.where(Elog[0] == e),np.where(theta[0] == th[t]))
    Sphi = Aproj[0][sel] # GRAND physical area projected along traj  [km2].

    ok_core = np.logical_and(np.logical_and(yc[0,sel]>-8500.,yc[0,sel]<8500.),np.logical_and(xc[0,sel]>-8000., xc[0,sel]<8000.))
    #ok_agg = np.logical_or(Nant_ew_agg[0][sel]>=nantmin,Nant_ns_agg[0][sel]>=nantmin)#,Nant_up_agg[0][sel]>=nantmin)
    #ok_cons = np.logical_or(Nant_ew_cons[0][sel]>=nantmin,Nant_ns_cons[0][sel]>=nantmin)#,Nant_up_cons[0][sel]>=nantmin)
    ok_agg = np.logical_and(np.logical_or(Nant_ew_agg[0][sel]>=nantmin,Nant_ns_agg[0][sel]>=nantmin),ok_core)#,Nant_up_agg[0][sel]>=nantmin)
    ok_cons = np.logical_and(np.logical_or(Nant_ew_cons[0][sel]>=nantmin,Nant_ns_cons[0][sel]>=nantmin),ok_core)#,Nant_up_cons[0][sel]>=nantmin)
    sel2 = np.where( (eny == np.round(10**e/1e18,2)) & (180-zen == th[t]) )  # In the other table
    Ragg = float(np.sum(ok_agg))/np.sum(Ntry[sel2]) #Ragg = float(np.sum(ok_agg))/np.size(ok_agg)
    Rcons = float(np.sum(ok_cons))/np.sum(Ntry[sel2]) #Rcons = float(np.sum(ok_cons))/np.size(ok_cons)

    Aeff_phi_agg[i,t] = Sphi[0]*Ragg   # Effective area integrated over phi
    Aeff_phi_cons[i,t] = Sphi[0]*Rcons   # Effective area integrated over phi
    Aeff_phi_opt[i,t] = Sphi[0]
    Aeff_phi_opt_200k[i,t] = Aproj_200k[0][sel][0]

	#pl.figure()
	#pl.show()

  aper_agg[i] = 2*pi*np.trapz(Aeff_phi_agg[i,:]*np.sin(th*pi/180),th*pi/180)
  aper_cons[i] = 2*pi*np.trapz(Aeff_phi_cons[i,:]*np.sin(th*pi/180),th*pi/180)
  aper_opt[i] = 2*pi*np.trapz(Aeff_phi_opt[i,:]*np.sin(th*pi/180),th*pi/180)

  expo_agg[i] = aper_agg[i]*dt_day
  expo_cons[i] = aper_cons[i]*dt_day
  expo_opt[i] = aper_opt[i]*dt_day

## Compute event rate
#day=24*3600

# 1e6 to transfer from km2 to m2
evt1d_agg = 1e6*(np.trapz(expo_agg[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_agg[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
evt1d_cons = 1e6*(np.trapz(expo_cons[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_cons[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
evt1d_opt = 1e6*(np.trapz(expo_opt[ind_le]*J1*pow(ee[ind_le],-gamma1),ee[ind_le])+np.trapz(expo_opt[ind_he]*J2*pow(ee[ind_he],-gamma2),ee[ind_he]))
evt1d_agg_diff = 1e9*1e6*np.concatenate([expo_agg[ind_le]*J1*pow(ee[ind_le],-gamma1),expo_agg[ind_he]*J2*pow(ee[ind_he],-gamma2)])  #GeV-1 day-1
evt1d_cons_diff = 1e9*1e6*np.concatenate([expo_cons[ind_le]*J1*pow(ee[ind_le],-gamma1),expo_cons[ind_he]*J2*pow(ee[ind_he],-gamma2)])

evt1d_agg_1e19 = 1e6*np.trapz(expo_agg[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])
evt1d_cons_1e19 = 1e6*np.trapz(expo_cons[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])
evt1d_opt_1e19 = 1e6*np.trapz(expo_opt[-2:]*J2*pow(ee[-2:],-gamma2),ee[-2:])

evt1d_agg_1e1819 = 1e6*np.trapz(expo_agg[1:3]*J2*pow(ee[1:3],-gamma2),ee[1:3])
evt1d_cons_1e1819 = 1e6*np.trapz(expo_cons[1:3]*J2*pow(ee[1:3],-gamma2),ee[1:3])
evt1d_opt_1e1819 = 1e6*np.trapz(expo_opt[1:3]*J2*pow(ee[1:3],-gamma2),ee[1:3])

evt1d_agg_lastbd = 1e6*np.trapz(np.array([expo_agg[-1],expo_agg[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))
evt1d_cons_lastbd = 1e6*np.trapz(np.array([expo_cons[-1],expo_cons[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))
evt1d_opt_lastbd = 1e6*np.trapz(np.array([expo_opt[-1],expo_opt[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))

#evt_eny_agg = expo_agg*1e6*J2*pow(ee,-gamma2) #/1e9
#evt_eny_cons = expo_cons*1e6*J2*pow(ee,-gamma2) #/1e9
#evt_eny_opt = expo_opt*1e6*J2*pow(ee,-gamma2) #/1e9
#evt1d_agg = np.trapz(evt_eny_agg,ee) #/1e9
#evt1d_cons = np.trapz(evt_eny_cons,ee) #/1e9
#evt1d_opt = np.trapz(evt_eny_opt,ee) #/1e9

print 'exposure_300 = ',expo_cons[-1]*np.arange(1,12)*365.*3600.*24./dt_day  *sec2yr

ndays_ref = 365.*1.  # duration of observation (days)
dt_ref = 3600.*24.*ndays_ref #in s
surf_ref = 200000 #10000 #in km2
rescale_surf = surf_ref/surf #Aeff_phi_opt_200k[0,0]/Aeff_phi_opt[0,0] #
rescale_dt = dt_ref/dt_day


print 1e6*np.trapz(expo_cons[0:2]*J2*pow(ee[2],-gamma2),ee[0:2])*5*rescale_surf *rescale_dt
print 1e6*np.trapz(expo_cons[2:4]*J2*pow(ee[2:4],-gamma2),ee[2:4])*5*rescale_surf *rescale_dt
print 1e6*np.trapz(np.array([expo_cons[-1],expo_cons[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))*5*rescale_surf *rescale_dt

Aeff_phi_agg = Aeff_phi_agg *rescale_surf
Aeff_phi_cons = Aeff_phi_cons *rescale_surf
Aeff_phi_opt = Aeff_phi_opt *rescale_surf
aper_agg = aper_agg *rescale_surf
aper_cons = aper_cons *rescale_surf
aper_opt = aper_opt *rescale_surf
expo_agg = expo_agg *rescale_surf *rescale_dt *sec2yr
expo_cons = expo_cons *rescale_surf *rescale_dt *sec2yr
expo_opt = expo_opt *rescale_surf *rescale_dt *sec2yr


print 'Energy:',ee,'eV'
#print 'Aperture (agressive):',aper_agg,'km2.sr   -   Exposure (aggressive):',expo_agg*sec2yr,'km2.sr.yr   -   Event rate (agressive):',evt_eny_agg
#print ' '
#print 'Aperture (conservative):',aper_cons,'km2.sr   -   Exposure (conservative):',expo_cons*sec2yr,'km2.sr.yr   -   Event rate (conservative):',evt_eny_cons
#print ' '
#print 'Optimal apperture:',aper_opt,'km2.sr   -   Exposure (optimal):',expo_opt*sec2yr,'km2.sr.yr   -   Event rate (optimal):',evt_eny_opt
#print ' '

print 'exposure_200k = ',expo_cons[-1]*np.arange(1,12)

print 'Expected daily event rate in 10^17-10^19.5eV (agressive):',evt1d_agg
print 'Expected daily event rate in 10^17-10^19.5eV (conservative):',evt1d_cons*rescale_surf
print 'Optimal daily event rate in 10^17-10^19.5eV:',evt1d_opt

print 'Expected daily event rate in 10^18-10^19eV (agressive):',evt1d_agg_1e1819
print 'Expected daily event rate in 10^18-10^19eV (conservative):',evt1d_cons_1e1819*rescale_surf
print 'Optimal daily event rate in 10^18-10^19eV:',evt1d_opt_1e1819

print 'Expected daily event rate for E>10^19eV (agressive):',evt1d_agg_1e19
print 'Expected daily event rate for E>10^19eV (conservative):',evt1d_cons_1e19*rescale_surf
print 'Optimal daily event rate for E>10^19eV:',evt1d_opt_1e19

print 'Expected daily event rate for E>10^19.5eV (agressive):',evt1d_agg_lastbd
print 'Expected daily event rate for E>10^19.5eV (conservative):',evt1d_cons_lastbd*rescale_surf
print 'Optimal daily event rate for E>10^19.5eV:',evt1d_opt_lastbd

evt1d_cons_diff = evt1d_cons_diff*rescale_surf


expo_auger = np.array([5400*dt_day]) #np.array([48000./(9.*365.*3600.*24.),48000./(9.*365.*3600.*24.)])
evt1d_cons_auger = 1e6*np.trapz(expo_auger*J2*pow(np.array([ee[-1:],1e20]),-gamma2),np.array([ee[-1:],1e20]))
#evt1d_cons_diff_auger = 1e9*1e6*expo_auger*J2*pow(ee[-1],-gamma2)

daily_event_rate_200k_1e19 = evt1d_cons_1e19 *rescale_surf #*rescale_dt#1e6*np.trapz(np.array([expo_cons[-1],expo_cons[-1]])*J2*pow(np.array([ee[-1],1e20]),-gamma2),np.array([ee[-1],1e20]))
print 'Expected daily event rate for E>10^19.5eV (conservative):',daily_event_rate_200k_1e19

print 'evt1d_cons_diff_auger:',evt1d_cons_auger
'''
for i in range(np.size(eu)):
  print '\n##E = 10^',eu[i],'eV'
  print Aeff_phi_agg[i,:]
  print Aeff_phi_cons[i,:]
  print Aeff_phi_opt[i,:]
'''
#print rescale_surf #*rescale_dt



exit()
#pl.figure(100)
#pl.semilogy(np.arange(1,11),expo_cons[-1]*np.arange(1,11)) #,'+-',lw=3,label='GRAND',mew=3, ms=10)
#print aper_agg,aper_cons
#print expo_agg*365.25*sec2yr,expo_cons*365.25*sec2yr
###############################################################
# Flux
pl.figure(17)
pl.loglog(leth,f1,lw=3,label='flux below ankle')
pl.loglog(heth,f2,lw=3,label='flux above ankle')
pl.xlabel('Energy (eV)')
pl.ylabel('Flux*E$^3$ (eV$^2$/s/sr/m$^2$)')
pl.legend(loc='best')
pl.title("TA flux - astro-ph 1511.07510")
pl.grid(True)

# Plot aperture, exposure
pl.figure(1)
#pl.loglog(ee,aper_agg,'+-',lw=3,label='Agressive',mew=3, ms=10)
#pl.loglog(ee,aper_cons,'+-',lw=3,label='Conservative',mew=3, ms=10)
#pl.loglog(ee,aper_opt,label='Optimal')

pl.loglog(ee,aper_cons,'+-',lw=3,label='GRAND',mew=3, ms=10)
pl.loglog([1e19,1e20],[5400,5400],'-',lw=3,label='Auger',mew=3,ms=10)
pl.xlabel('Energy (eV)')
pl.ylabel('Aperture (km$^2$.sr)')
pl.grid(True)
pl.legend(loc='best')
#pl.title("GRANDproto300 CR aperture")
pl.title("GRAND CR aperture")

pl.figure(2)
#pl.loglog(ee,expo_agg,'+-',lw=3,label='Agressive',mew=3, ms=10) #*365.25
#pl.loglog(ee,expo_cons,'+-',lw=3,label='Conservative',mew=3, ms=10) #*365.25
#pl.loglog(ee,expo_opt*365.25,label='Optimal')

pl.loglog(ee,expo_cons,'+-',lw=3,label='GRAND',mew=3, ms=10) #*365.25
pl.loglog([1e19,1e20],[48000,48000],'-',lw=3,label='Auger',mew=3,ms=10)
pl.xlabel('Energy (eV)')
pl.ylabel('Exposure (km$^2$.sr.s)')
pl.grid(True)
pl.legend(loc='best')
#pl.title("GRANDProto300 CR exposure")
pl.title("GRAND CR exposure")

pl.figure(3)
pl.semilogy(th,Aeff_phi_opt[0,:],lw=3,label='Optimal')
pl.semilogy(th,Aeff_phi_opt_200k[0,:],lw=3,label='Optimal')
for ie in range(np.size(ee)):
	pl.semilogy(th,Aeff_phi_agg[ie,:],'+-',lw=3,label='Agressive - '+str(np.round(ee[ie]/1e18,1))+' EeV',mew=3, ms=10)
pl.xlabel('Theta (deg)')
pl.ylabel('$<A_{eff}>_\phi$ (km$^2$)')
#pl.ylim([2,7000])
pl.grid(True)
pl.legend(loc='lower right',fontsize=11)
pl.title("GRANDProto300 CR effective area (agressive case)")

pl.figure(4)
pl.semilogy(th,Aeff_phi_opt[0,:],lw=3,label='Optimal')
pl.semilogy(th,Aeff_phi_opt_200k[0,:],lw=3,label='Optimal')
for ie in range(np.size(ee)):
	pl.semilogy(th,Aeff_phi_cons[ie,:],'+-',lw=3,label='Conservative - '+str(np.round(ee[ie]/1e18,1))+' EeV',mew=3, ms=10)
pl.xlabel('Theta (deg)')
pl.ylabel('$<A_{eff}>_\phi$ (km$^2$)')
#pl.ylim([2,7000])
pl.grid(True)
pl.legend(loc='lower right',fontsize=11)
pl.title("GRANDproto300 CR effective area (conservative case)")

pl.figure(5)
#pl.loglog(ee,evt1d_agg_diff,'+-',lw=3,label='Agressive',mew=3, ms=10)
pl.loglog(ee,evt1d_cons_diff,'+-',lw=3,label='Conservative',mew=3, ms=10)
#pl.loglog(ee,evt_eny_opt*sec2yr,label='Optimal')
pl.xlabel('Energy (eV)')
pl.ylabel('dN/dE/dt (GeV$^{-1}$ day$^{-1}$)')
pl.grid(True)
#pl.legend(loc='best')
pl.title("GRANDProto300 CR event rate")
pl.show()
