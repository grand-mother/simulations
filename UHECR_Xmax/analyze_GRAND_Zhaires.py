### This script is friendly provided by Stijn Buitink!
### It was already adapted to perform Xmax reco for SKA-low.
### Now adapted to perform Xmax reco for GRAND.

import numpy as np
from Tkinter import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import cPickle
import scipy.interpolate as intp
import scipy.optimize as opt
import random
import NoisePower
from scipy.optimize import least_squares

import Remove_outliers


# Parabola model, b=sx, c=sy
def model(x, a, b, c):# model(x, a=1.0, b=700., c=1.65e9)
  return a*(x-b)**2 +c

def GetPower_UVW(amp, power, cx,cy,zen,az):
   #first transform to vector
   vec=np.sign(amp)*np.sqrt(power)
   #no get UVW coordinates
   vec2=GetUVW(vec,cx,cy,zen,az)
   #transform to power
   return np.square(vec2)

def GetXYZ(pos, zen, az):
   inc=152.95*np.pi/180. #no raspass used #(180.-27.05)*np.pi #-60.0/180.*np.pi #magnetic field direction on GRAND site
   B = np.array([0,np.cos(inc),-np.sin(inc)])
   v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
   vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
   return pos[0]*vxB+pos[1]*vxvxB+pos[2]*v

def GetUVW(pos, cx, cy,cz, zen, az):
   
   relpos = pos-np.array([cx,cy,cz])
   inc=152.95*np.pi/180.
   phigeo = 0.*np.pi/180.

   B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) #from oliviers script including phigeo
   B=B/np.linalg.norm(B)
   v = np.array([np.cos(az)*np.sin(zen),np.sin(az)*np.sin(zen),np.cos(zen)]) # or *-1: change the direction
   v=v/np.linalg.norm(v)
   vxB = np.cross(v,B) #np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.cross(v,vxB) #np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
   vxvxB = vxvxB/np.linalg.norm(vxvxB)
   return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T





def reverseAnalysis(simfile, eventno, simevent, eventno2,  flagging=True, plots=True, 
                    verbose=True, outfile="reverse", simmode=True, saveplt=False, 
                    showplt=True,  outputfolder=True, SKAmap=True):
  
   print "\nStart reverse analysis"
   
   realxmax=0
   random.seed()
   
   g=open(simfile,'r')
   siminfo = cPickle.load(g)
   g.close()

   zenith=siminfo['zenith']
   azimuth=siminfo['azimuth']
   hillas=siminfo['hillas'] # parallel version
   sim_antenna_position=siminfo['antenna_position'] # just x and y??
   primary=siminfo['primary']
   lowco=siminfo['lowco'][0]
   hico=siminfo['hico'][0]
#   sim_power=siminfo['filteredpower']
   sim_power=siminfo['power41']

   del siminfo

   # Simulated star shape has to be centralised (to be included-done). Actual array can than be shifted by any number
   # Here you define core position of simulation in meters wrt to array core NOTE correct that
   core_x=0 #np.mean(sim_antenna_position[0][:,0]) #50
   core_y=0 #np.mean(sim_antenna_position[0][:,1]) #100  # - = to the right0 #0#
   core_z=0 #np.mean(sim_antenna_position[0][:,2]) #
   print "core off by ", core_x, core_y, core_z

   rand_core_x = random.gauss(0,500)
   rand_core_y = random.gauss(0,500)

   # Add an offset to correct for the shower core in the simulations
   corr_x=np.mean(sim_antenna_position[0][:,0])+rand_core_x #50
   corr_y=np.mean(sim_antenna_position[0][:,1])+rand_core_y #100  # - = to the right0 #0#
   corr_z=np.mean(sim_antenna_position[0][:,2]) #
   print "corrected by ", corr_x, corr_y, corr_z
  
   sim_tot_power=np.sum(sim_power,axis=2)
   
   nsim=sim_tot_power.shape[0] # rows = number of Simulations
   nsimant=sim_tot_power.shape[1] #columns = number of Antennas
   
   for i in np.arange(int(nsim)):
    sim_antenna_position[i][:,0]= sim_antenna_position[i][:,0] -corr_x
    sim_antenna_position[i][:,1]= sim_antenna_position[i][:,1] -corr_y
    sim_antenna_position[i][:,2]= sim_antenna_position[i][:,2] -corr_z
   
   data_zenith=zenith[0]
   data_azimuth=azimuth[0]

   pos_sim_UVW = np.zeros([nsim,nsimant,3])
   
   #### Interpolation on simulated antenna grid starts
   # Interpolations based on RBF functions:
   
   rbf= np.ndarray([nsim],dtype=object)

   def func_call(self,r):
       return np.log((r/self.epsilon)**2+1)

   for i in np.arange(int(nsim)):
      pos_sim_UVW[i,:,:] = GetUVW(sim_antenna_position[i][:,:],0,0,0, zenith[0], azimuth[0]) #transform to shower plane
      
      selection=np.array(np.isfinite(sim_power[i][:,0]))*np.array(np.isfinite(sim_power[i][:,1]))
      
      rbf[i] = intp.Rbf(pos_sim_UVW[i,selection,0], pos_sim_UVW[i,selection,1], sim_tot_power[i,selection],smooth=0,function=func_call)

   # 2D radio fit function
   def radio_fitfunc(f,dpos,n,cx,cy,az,zen):
       pos_ant_UVW = GetUVW(dpos, cx, cy,0, zen, az)
       interp_power = rbf[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
       bool1=interp_power>=0
       interp_power *= bool1
       return f*interp_power

   # radio fit p=[pratio,xoff,yoff]
   radio_errfunc = lambda p,lofar_pow,lofar_err, lofar_pow01,lofar_err01, lofar_pos,lora_data, lora_err, lora_pos, cx, cy, az, zen, n: (radio_fitfunc(p[0],lofar_pos,n,cx+p[1],cy+p[2],az,zen) - lofar_pow)/lofar_err

   # combined error function p = [pratio,dratio,xoff,yoff]
   def combined_errfunc(p,*args):
     a = radio_errfunc([p[0],p[2],p[3]],*args)
     b = radio_errfunc([p[0],p[2],p[3]],*args)*0. #lora_errfunc_restrained([p[1],p[2],p[3]],*args)
     c = np.concatenate((a,b))
     return c
 

 
 
   ##################################################################
 
 
 
   ## Now create the "fake" data:
   realxmax=hillas[simevent] # parallel version
   print "realxmax ", realxmax
   
   ## read in GRAND positions!
   skafile=open(SKAmap,'r')
   lines=skafile.readlines()


   nant=len(lines) # number of GRAND antennas
   ska_x=np.zeros([nant])
   ska_y=np.zeros([nant])
   ska_z=np.zeros([nant])
   for i in np.arange(nant): # - corr to be deleted if array around (000)
       ska_x[i]=lines[i].split()[0]
       ska_y[i]=lines[i].split()[1]
       ska_z[i]=lines[i].split()[2] 
       ska_x[i]= ska_x[i]- corr_x
       ska_y[i]= ska_y[i]- corr_y
       ska_z[i]= ska_z[i]- corr_z
   skafile.close()
   del skafile


   ### ANTENNA selection: ex. within 300m around core and just every 8 for SKA
   # Play around with that !
   # NOTE: set radius in order to include all antennas

   ant_sel=np.ones([nant],dtype=bool)

   ant_sel=(np.sqrt(np.square(ska_x-core_x)+np.square(ska_y-core_y)+np.square(ska_z-core_z))<40000)*(np.arange(0,nant)%1==0)
#   print np.sqrt(np.square(ska_x-core_x)+np.square(ska_y-core_y)+np.square(ska_z-core_z)), (np.arange(0,nant)%2==0)

   nsel_ant = np.sum(ant_sel)
   #print "selected: ", nsel_ant # number of antennas used in the analysis
   positions=np.zeros([nant,3])
   positions[:,0]=ska_x
   positions[:,1]=ska_y
   positions[:,2]=ska_z #np.ones([nant])*370.
   del ska_x, ska_y, ska_z
   pos_ant_UVW = GetUVW(positions, core_x, core_y, core_z,data_zenith, data_azimuth) # rotates antennas in vxB-vxvxB
   
   
   
   
   ### Interpolated footprint used to set the power to each GRAND antenna position
   dtotpower=np.zeros([len(pos_ant_UVW[:,0])])
   for i in np.arange(len(pos_ant_UVW[:,0])):
      dtotpower[i]=rbf[simevent](pos_ant_UVW[i,0],pos_ant_UVW[i,1]) # give a power value to all GRAND antenna position based on all simulation
      if dtotpower[i]<0:
           dtotpower[i]=0

   
   ### NOISE section: has to be deleted for a more realistic one, ATTENTION units not correct

   noisefraction=0.01
   print "Noiselevel ", noisefraction
   print "galactic noise 0.01"

   # always a specific fraction of the power as noise for the every single antenna + 1% of the power of the antenna with the highest signal for galactic noise
   #dsigmatot=np.full((len(pos_ant_UVW[:,0])),1.0**2)

#   dsigmatot=dtotpower*noisefraction+0.01*np.max(dtotpower) #change this to a more realistic!!!

#   print "Pnoise", Pnoise

   nsample = 512. # number of bins for intergration along time traces
   Pnoise = NoisePower.Calc_Noise(zenith[0],azimuth[0],lowco,hico)*1.
   dsigmatot=np.sqrt(4.*Pnoise*dtotpower + 2.*nsample*Pnoise**2) # New noise model, for voltage traces

   dpower=np.ones([nant,2])
   dpower[:,0]=dtotpower/2
   dpower[:,1]=dtotpower/2
   dsigma=np.ones([nant,2])
   dsigma[:,0]=dsigmatot/2
   dsigma[:,1]=dsigmatot/2
   simpower=np.zeros([nant,nsim])



   for i in np.arange(nant):
      dtotpower[i]=dtotpower[i]+random.gauss(0,dsigmatot[i]) # Add the noise to the power

   # Threshold value for antenna triggering (power, rough estimate)
#   thresh = 150.**2
#      if dtotpower[i]<thresh:
#            dtotpower[i]=0

   print "Mean noise, Max noise, Max power (muV/m)^2", np.mean(dsigmatot[ant_sel]),np.max(dsigmatot), np.max(dtotpower)
   

   ################ PARTICLE DETECTOR #insert particle detector coordinates here
   nstations=1 
   lora_positions=np.zeros([nstations,3])
   lora_positions[:,0]=10
   lora_positions[:,1]=10

   lora_dens = np.ones([nstations,1])*0.#lora_fitfunc(1.0,lora_positions,0, 0,data_azimuth,data_zenith,simevent)
   lora_err = np.sqrt(lora_dens)

   
   ##### STARTING CHI2 procedure

   combchi2=np.zeros([nsim])
   radiochi2=np.zeros([nsim])
   lorachi2=np.zeros([nsim])
   p_ratio=np.zeros([nsim])
   d_ratio=np.zeros([nsim])
   xoffset=np.zeros([nsim])
   yoffset=np.zeros([nsim])
   
   
   # Getting Chi2
   def FitRoutine(fit_args,iterations=1, dosim=False):
      paramset=np.ndarray([iterations,4])
      itchi2=np.zeros(iterations)
      for j in np.arange(iterations):
         initguess=[1,1,200*np.random.rand()-100,200*np.random.rand()-100] #CHANGE VALUES ACCRODING TO MODE!!!
         paramset[j], covar = opt.leastsq(combined_errfunc,initguess, args=fit_args)
         itchi2[j]=np.sum(np.square(combined_errfunc(paramset[j],*fit_args)))
         #if (verbose): print itchi2[j], paramset[j]
      bestiter=np.argmin(itchi2)
      fitparam=paramset[bestiter]
      chi2_comb=np.sum(np.square(combined_errfunc(fitparam,*fit_args)))
      chi2_rad=np.sum(np.square(radio_errfunc([fitparam[0],fitparam[2],fitparam[3]],*fit_args)))
      chi2_part=0.#np.sum(np.square(lora_errfunc_restrained([fitparam[1],fitparam[2],fitparam[3]],*fit_args)))
      return chi2_comb, chi2_rad, chi2_part, fitparam[0], fitparam[1], fitparam[2], fitparam[3]
      del chi2_comb,chi2_rad , chi2_part, paramset, itchi2
   
   print "len dtotpower",len(dtotpower)

   niterations=1
   for i in np.arange(nsim):
      fit_args=(dtotpower[ant_sel],dsigmatot[ant_sel],dpower[ant_sel],dsigma[ant_sel],positions[ant_sel],lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,i)
      combchi2[i], radiochi2[i], lorachi2[i], p_ratio[i], d_ratio[i], xoffset[i], yoffset[i]= FitRoutine(fit_args,niterations, simmode)

   nsel_ant = np.sum(ant_sel)
   print "Antennas flagged: ", nsel_ant, " out of ", nant

   ndf_comb = nsel_ant+nstations-5
   ndf_radio = nsel_ant-4


   if (simmode): 
      combchi2[simevent]=1000*ndf_comb #NOTE: penalty for actual event used to make fake data
      radiochi2[simevent]=1000*ndf_radio
   bestsim=np.argmin(combchi2)
   xoff=xoffset[bestsim]
   yoff=yoffset[bestsim]
   
   #ADDITIONAL FLAGGING
   fit_args=(dtotpower,dsigmatot,dpower,dsigma,positions,lora_dens,lora_err,
             lora_positions,core_x,core_y,data_azimuth,data_zenith,bestsim)
  

   # Second round of fitting after flagging of antennas:
   for i in np.arange(nsim):
      fit_args=(dtotpower[ant_sel],dsigmatot[ant_sel],dpower[ant_sel],dsigma[ant_sel],positions[ant_sel],lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,i)
      combchi2[i], radiochi2[i], lorachi2[i], p_ratio[i], d_ratio[i], xoffset[i], yoffset[i]= FitRoutine(fit_args,niterations, simmode)
      simpower[:,i]=radio_fitfunc(p_ratio[i],positions,i,core_x+xoffset[i],core_y+yoffset[i],azimuth[0],zenith[0])
   
   pos_ant_UVW = GetUVW(positions, core_x+xoff, core_y+yoff,core_z, data_zenith, data_azimuth)
   axdist_ant = np.sqrt(pos_ant_UVW[:,0]*pos_ant_UVW[:,0]+pos_ant_UVW[:,1]*pos_ant_UVW[:,1])

   
   print "Starting fit procedure"
   
   # Select points for the fit:
   
   fit_selection=np.ones(nsim) # all points that have lower chi2 values on one side only

#   fit_selection=np.zeros(nsim) # all points that have lower chi2 values on one side only
#   for i in np.arange(nsim):
#      if (np.sum(combchi2[int(hillas[:]>hillas[i])]<combchi2[i])==0): fit_selection[i]=fit_selection[i]+1 # parallel
#      if (np.sum(combchi2[int((hillas[:]<hillas[i]))]<combchi2[i])==0): fit_selection[i]=fit_selection[i]+1 # parallel
   
   if (simmode): fit_selection[simevent]=0 # for fake data

   
   # Proton
   h_p=[]
   c_p=[]
   for i in np.arange(nsim):
     if (primary[i] <20):
       h_p.append(hillas[i])# for plotting
       c_p.append(combchi2[i])# for plotting

   # Iron
   h_Fe=[]
   c_Fe=[]
   for i in np.arange(nsim):
     if (primary[i] >20):
        h_Fe.append(hillas[i]) # for plotting
        c_Fe.append(combchi2[i])# for plotting
   
   # All: these simulation will go into the fit procedure
   h=[]
   c=[]
   for i in np.arange(nsim):
       if(fit_selection[i]>0):
          h.append(hillas[i])
          c.append(combchi2[i])

   #######################################
   ### Fitting selection and procedure ###
   #######################################


   t_new = np.linspace(min(h[:]), max(h[:]), 50)
   def get_flipped(y_data, y_model):
      flipped = y_model - y_data
      flipped[flipped > 0] = 0
      return flipped

   def flipped_resid(pars, x, y):
    #For every iteration, everything above the currently proposed
    #curve is going to be mirrored down, so that the next iterations
    #is going to progressively shift downwards.
      y_model = model(x, *pars)
      flipped = get_flipped(y, y_model)
      resid = np.square(y + flipped - y_model)
      return np.nan_to_num(resid)

   guesses = [2.0e-4, 800., 10.0]
   
   res1 = least_squares(flipped_resid, guesses, args = (h, c/(ndf_comb+1e-25)), 
                        jac='2-point', bounds=(0., np.inf), method='trf')

   fit_pars = res1.x
   y_fit = model(t_new, *fit_pars)
   Xreco = fit_pars[1]
   
   print 'Edge of point cloud: Curvature, Xmax_min, Chi2_min ', fit_pars


   m1, m2 = np.array(h), np.array(c/(ndf_comb+1e-25))
   xmin = m1.min()-10.
   xmax = m1.max()+10.

   Xcol = np.arange(xmin,xmax+(xmax-xmin)/9.,(xmax-xmin)/9.)                

   # Fit selection: as a function of point distribution in each Xmax bin
   hfit_new = np.array([])
   cfit_new = np.array([])
   for i in range(len(m1)):
       ni = int((m1[i]-xmin)/(Xcol[1]-Xcol[0]))
       minChi2 = np.min( m2[ (m1>=Xcol[ni]) & (m1<Xcol[ni+1]) ] )
       stdChi2 = np.std( m2[ (m1>=Xcol[ni]) & (m1<Xcol[ni+1]) ] )
       if (m2[i] <= minChi2+stdChi2):
           hfit_new=np.append(hfit_new, m1[i])
           cfit_new=np.append(cfit_new, m2[i])

   print('Simulations selected for fitting: %i' %len(hfit_new) + '/%i' %len(hillas))
   
   try:
      popt, pcov = opt.curve_fit(model, hfit_new, cfit_new,
                    bounds=([0.,0.,0.],[np.inf,np.inf,np.inf])) # with bounds
   
      print 'Curve Fit with selection: Curvature, Xmax_min, Chi2_min ', popt

      Xreco = popt[1]
      y_fit = model(t_new, *popt)

   except (RuntimeError, TypeError, NameError, UnboundLocalError) as exc:
       print(exc)
       Xreco=hillas[bestsim]

   

#==============================================================================
# Figures
#==============================================================================

   ### Figure with interpolated footprint in vxB and simulated star shape with signal

   fig=plt.figure(1,figsize=(8,8)) 

   dist_scale = np.max(abs(pos_sim_UVW[bestsim,:,0]))
   ti = np.linspace(-dist_scale, dist_scale, 150)
   XI, YI = np.meshgrid(ti, ti)
   ZI = rbf[bestsim](XI, YI)*p_ratio[bestsim] # gives you a plotting grid for plotting the interpolated footprint
   maxp = np.max([np.max(dtotpower),np.max(ZI)]) # if =1, than there is no scaling of the signal distribution to the highest power1 #
   
   plt.imshow(ZI/maxp,vmax=1, vmin=-0.03, interpolation='bilinear', extent=(-dist_scale,dist_scale,-dist_scale,dist_scale), origin='lower left',cmap=cm.gnuplot2_r) # interpolation
   im_vxB=plt.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],15,sim_tot_power[bestsim,:]*p_ratio[bestsim]/maxp, vmax=1, vmin=-0.03,cmap=cm.gnuplot2_r, lw = 0.1 ) # starshape pattern
   
   fig.colorbar(im_vxB) # im, ax=ax
   plt.xlabel(r"Position along ${\bf v} \times \,{\bf B}$ axis (m)", fontsize=16)
   plt.ylabel(r"Position along ${\bf v} \times \, ({\bf v} \times \,{\bf B})$ axis (m)", fontsize=16)
   name=outputfolder+"GRAND_eventview_{0}_{1}.png".format(eventno, simevent)
   plt.savefig(name)

   plt.close()
   
   
   ### Figure with LDF

   fig=plt.figure(2) 

   maxp = np.max([np.max(dtotpower),np.max(simpower[:,bestsim])]) #1 #
   dataplt=plt.errorbar(axdist_ant[ant_sel],dtotpower[ant_sel]/maxp,dsigmatot[ant_sel]/maxp,linestyle='',marker='o',color='r')
   simplt=plt.plot(axdist_ant[ant_sel],simpower[ant_sel,bestsim]/maxp,linestyle='',marker='s',color='b')
   plt.ylabel("Total power (a.u.)", fontsize=16)
   plt.xlabel("Distance (m)", fontsize=16)
   plt.ylim((0,max(simpower[ant_sel,bestsim]/maxp)+0.2))
   plt.xlim(0,4000)
   plt.legend([dataplt[0],simplt[0]],["Fake data","ZHAireS simulations"], numpoints=1, frameon=False)
   name=outputfolder+"GRAND_ldf_{0}_{1}.png".format(eventno, simevent)
   plt.savefig(name)
 
   plt.close()


   ### Chi2 Xmax and fitted parabola
   
   fig=plt.figure(4,figsize=(8,8))  
   
   plt.ylim(0,2)

   plt.ylabel(r"$\chi^2$/ ndf", fontsize=16)
   plt.xlabel(r"$X_{\mathrm{max}}$ (g/cm$^2$)", fontsize=16)
 
   plt.xlim(500,1100)

   plt.plot(t_new, y_fit)

   plt.scatter(h_p,c_p[:]/(ndf_comb+1e-25),50,color='r')
   plt.scatter(h_Fe,c_Fe[:]/(ndf_comb+1e-25),50,color='b', marker='s')

   plt.scatter(hfit_new, cfit_new,50, color='#2ca02c', marker='x')

   plt.legend(["Fitted Curve","ZHAireS Simulations (proton)","ZHAireS Simulations (iron)"],
              fontsize=14,frameon=False)
   name=outputfolder+"GRAND_xmaxcurve_{0}_ev{1}_realXmax{2}.png".format(eventno, simevent, realxmax)
   plt.savefig(name)  
   plt.close()

   
   
   ### GRAND eventview on ground
   # Figure with GRAND antenna positions and their power on ground

   fig=plt.figure(7,figsize=(8,6)) 
   plt.plot(core_x, core_y, '+', ms=6)
   dist_scale3 = 400 
   ti3 = np.linspace(-dist_scale3, dist_scale3, 150)
   XI3, YI3 = np.meshgrid(ti3, ti3)

   ZI3 = rbf[bestsim](XI3, YI3)*p_ratio[bestsim] # For the interpolated footprint, but not needed here, at some point you need it to get the right color bar scale
   maxp = np.max([np.max(dtotpower),np.max(ZI3)])
   cmap = matplotlib.cm.get_cmap("jet")

   im2=plt.scatter(positions[:,0],positions[:,1],10,dtotpower/maxp,
                   vmax=1, vmin=-0.03,cmap=cmap, lw = 0 )

   dist_scale=np.max([abs(positions[:,0]), abs(positions[:,1])])
   plt.xlim((-dist_scale,dist_scale))
   plt.ylim((-dist_scale,dist_scale))
   plt.xlabel("West-East (m)", fontsize=16)
   plt.ylabel("South-North (m)", fontsize=16)
   fig.colorbar( im2) # im, ax=ax

   name=outputfolder+"GRAND_eventviewground_{0}_{1}.png".format(eventno, simevent)
   plt.savefig(name)
   plt.close()


   ### Starshape on ground for comparison
   
   fig=plt.figure(5,figsize=(8,8))
   ti = np.linspace(-dist_scale, dist_scale, 150)
   XI, YI = np.meshgrid(ti, ti)
   selection = np.unique(np.array(np.isfinite(sim_tot_power[bestsim,:]))*np.arange(nsimant))
   ZI=np.zeros([150,150])
   for ii in np.arange(150):
      for jj in np.arange(150):
         gridpos = GetUVW(dist_scale*np.array([2*(ii/150.)-1,2*(jj/150.)-1,0]),0,0,0, data_zenith,data_azimuth)
         if (gridpos[0]*gridpos[0]+gridpos[1]*gridpos[1]<500*500): ZI[ii,jj] = rbf[bestsim](gridpos[0],gridpos[1])*p_ratio[bestsim]
   maxp = np.max([np.max(dtotpower),np.max(ZI)])
      
   plt.scatter(sim_antenna_position[bestsim][:,0]+core_x,sim_antenna_position[bestsim][:,1]+core_y,10, c=sim_tot_power[bestsim,:]*p_ratio[bestsim],vmax=maxp, vmin=0,cmap=cm.gnuplot2_r)

   plt.xlim((-dist_scale,dist_scale))
   plt.ylim((-dist_scale,dist_scale))
   plt.xlabel("West-East (m)")
   plt.ylabel("South-North (m)")
   name=outputfolder+"GRAND_ground_{0}.png".format(eventno)
   plt.savefig(name)  
 

   if Xreco<0.:
       print "\n"
       print "WARNING: fit not working correctly"
       print "\n"
               
   xmaxreco=0
   return Xreco, realxmax, hillas[bestsim], xmaxreco,c, h, (ndf_comb+1e-25), primary[simevent]


