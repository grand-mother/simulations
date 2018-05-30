#import matplotlib
#matplotlib.use('Agg')

import os, glob
import pylab as pl
#import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
import sys
import numpy as np
import linecache
import StringIO
#from scipy.optimize import curve_fit
#from scipy.misc import factorial
import scipy.stats as stats
pl.rc('text', usetex=True)
pl.rc('font', family='serif')

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<0 or len(sys.argv)>2):
    print("""\
        This script produces the histogram comparing the simulations performed with RadioMorphing and ZHAireS based on Efield files.

        Usage:  python histo_RM.py [1,2,3 or 4 depending on the wanted histograms]
        """)
    sys.exit(1)

prefix = 'a'
figdir = '/Users/nrenault/Desktop/GRAND/RadioMorphing/fig/'
DISPLAY = True

############################################################################################################
############################################################################################################
############################################################################################################
def read_showers(showers):
    zen_sh = []
    az_sh = []
    eny_sh = []
    injh_sh = []
    showerID_sh = []
    for shower in showers:
        showerID_sh.append(shower.replace('/','.').split('.')[-2])
        input_file_path = glob.glob(shower+'/*.inp')[0]
        zen,azim,energy,injh,primarytype = inputfromtxt(input_file_path)
        zen_sh.append(zen)
        az_sh.append(azim)
        eny_sh.append(energy)
        injh_sh.append(injh)

    showerID_sh = np.array(showerID_sh)
    zen_sh = np.array(zen_sh)
    az_sh = np.array(az_sh)
    eny_sh = np.array(eny_sh)
    injh_sh = np.array(injh_sh)

    return showerID_sh,zen_sh,az_sh,eny_sh,injh_sh

############################################################################################################
def find_ind(showerID_sh,zen_sh,az_sh,eny_sh,injh_sh,zen_bin,az_bin,eny_bin,injh_bin):
    #digitize provide the index of the upper boundary of the bin each value belongs to
    zen_ind = np.digitize(zen_sh, zen_bin)
    az_ind = np.digitize(az_sh, az_bin)
    eny_ind = np.digitize(eny_sh, eny_bin)
    injh_ind = np.digitize(injh_sh, injh_bin)

    return zen_ind,az_ind,eny_ind,injh_ind

############################################################################################################
def build_binning():
    zen_bin = np.arange(84.,90.*1.000001,1.)
    az_bin = np.arange(0.,360.*1.000001,22.5)
    eny_bin = 10**np.arange(np.log10(5e16),np.log10(1e19)*1.000001,0.5)
    injh_bin = np.arange(1500.,6000.*1.000001,250.)

    return zen_bin,az_bin,eny_bin,injh_bin

############################################################################################################
def inputfromtxt(input_file_path):
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

############################################################################################################
def Efield_zhaires(wkdir):
    try:
      fname = wkdir+'/split/antpos.dat' #ZHAires_output
      a = np.loadtxt(fname, dtype='float', comments='#')
    except:
      fname = wkdir+'/antpos.dat' #ZHAires_output
      a = np.loadtxt(fname, dtype='float', comments='#')

    if not(os.path.exists(fname)):
      print("no antpos.dat file")
      print("No antenna within the shower footprint")
      sys.exit(1)

    ###
    # First load file
    x0_zhaires = a[:,0] #X = S->N
    y0_zhaires = a[:,1] #Y = E->W
    z0_zhaires = a[:,2]
    Nant_zhaires = np.size(x0_zhaires)

    Ampx_zhaires = np.zeros((Nant_zhaires,1))
    Ampy_zhaires = np.zeros((Nant_zhaires,1))
    Ampz_zhaires = np.zeros((Nant_zhaires,1))
    for iant in range(0,Nant_zhaires):
        try:
          filename = wkdir+"/split/"+prefix+str(iant)+".trace"
          if not(os.path.exists(filename)):
            filename = wkdir+"/"+prefix+str(iant)+".trace"
          b = np.loadtxt(filename, dtype='float', comments='#')
          ti_zhaires = b[:,0]
          Ampx_zhaires[iant] = max(abs(b[:,1]))
          Ampy_zhaires[iant] = max(abs(b[:,2]))
          Ampz_zhaires[iant] = max(abs(b[:,3]))
        except:
          pass

    # Get amplitude min/max
    maxEx_zhaires = np.amax(Ampx_zhaires)
    maxEy_zhaires = np.amax(Ampy_zhaires)
    maxEz_zhaires = np.amax(Ampz_zhaires)
    minEx_zhaires = np.amin(Ampx_zhaires)
    minEy_zhaires = np.amin(Ampy_zhaires)
    minEz_zhaires = np.amin(Ampz_zhaires)

    return Ampx_zhaires,Ampy_zhaires,Ampz_zhaires,x0_zhaires,y0_zhaires,z0_zhaires

############################################################################################################
def Efield_RM(wkdir):
    try:
      fname = wkdir+'/split/antpos.dat' #ZHAires_output
      a = np.loadtxt(fname, dtype='float', comments='#')
    except:
      fname = wkdir+'/antpos.dat' #ZHAires_output
      a = np.loadtxt(fname, dtype='float', comments='#')

    if not(os.path.exists(fname)):
      print("no antpos.dat file")
      print("No antenna within the shower footprint")
      sys.exit(1)

    ###
    # First load file
    x0_RM = a[:,0] #X = S->N
    y0_RM = a[:,1] #Y = E->W
    z0_RM = a[:,2]
    Nant_RM = np.size(x0_RM)

    Ampx_RM = np.zeros((Nant_RM,1))
    Ampy_RM = np.zeros((Nant_RM,1))
    Ampz_RM = np.zeros((Nant_RM,1))
    for iant in range(0,Nant_RM):
        try:
          filename = wkdir+"/"+prefix+str(iant)+".trace"
          if not(os.path.exists(filename)):
            filename = wkdir+"/"+prefix+str(iant)+".trace"
          b = np.loadtxt(filename, dtype='float', comments='#')
          ti_RM = b[:,0]
          Ampx_RM[iant] = max(abs(b[:,1]))
          Ampy_RM[iant] = max(abs(b[:,2]))
          Ampz_RM[iant] = max(abs(b[:,3]))
        except:
          pass

    # Get amplitude min/max
    maxEx_RM = np.amax(Ampx_RM)
    maxEy_RM = np.amax(Ampy_RM)
    maxEz_RM = np.amax(Ampz_RM)
    minEx_RM = np.amin(Ampx_RM)
    minEy_RM = np.amin(Ampy_RM)
    minEz_RM = np.amin(Ampz_RM)

    return Ampx_RM,Ampy_RM,Ampz_RM,x0_RM,y0_RM,z0_RM

############################################################################################################
def histo_DEpp(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires):

    DAmpx_rel_all = []
    DAmpy_rel_all = []
    DAmpz_rel_all = []
    DAmpx_all = []
    DAmpy_all = []
    DAmpz_all = []
    for iserie in Nserie:
        for ie in range(len(Etarget)):
            ratioE = float(Etarget[ie])/1E+19

            for iDd in Dd:
                for ialfa in range(np.size(alfa)):
                    if alfa[ialfa]>0 or (alfa[ialfa]==0 and iDd==20000):
                        folder_zhaires=ROOT_zhaires+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+'/'
                        folder_RM=ROOT_RM+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'/'
                        showers = glob.glob(folder_zhaires+'[0-9][0-9]*[0-9]/')
                        #showers=showers[0:11]

                        #Start loop on shower
                        for shower in showers:
                            showerID = shower.replace('/','.').split('.')[-2]
                            trace_files=glob.glob(shower+'/split/*.trace')
                            Nant=len(trace_files)
                            if Nant>0:
                                Ampx_RM,Ampy_RM,Ampz_RM,x0_RM,y0_RM,z0_RM = Efield_RM(folder_RM+showerID+'/')
                                Ampx_zhaires,Ampy_zhaires,Ampz_zhaires,x0_zhaires,y0_zhaires,z0_zhaires = Efield_zhaires(shower)

                                ind = np.where(Ampy_RM==0)
                                Ampx_RM = np.delete(Ampx_RM,ind)
                                Ampy_RM = np.delete(Ampy_RM,ind)
                                Ampz_RM = np.delete(Ampz_RM,ind)
                                Ampx_zhaires = np.delete(Ampx_zhaires,ind)
                                Ampy_zhaires = np.delete(Ampy_zhaires,ind)
                                Ampz_zhaires = np.delete(Ampz_zhaires,ind)

                                ## Computation part
                                DAmpx = Ampx_RM-Ampx_zhaires
                                DAmpy = Ampy_RM-Ampy_zhaires
                                DAmpz = Ampz_RM-Ampz_zhaires
                                DAmpx_rel = (Ampx_RM-Ampx_zhaires)/Ampx_zhaires
                                DAmpy_rel = (Ampy_RM-Ampy_zhaires)/Ampy_zhaires
                                DAmpz_rel = (Ampz_RM-Ampz_zhaires)/Ampz_zhaires

                                DAmpx_rel = np.reshape(DAmpx_rel,np.size(DAmpx_rel))
                                DAmpy_rel = np.reshape(DAmpy_rel,np.size(DAmpy_rel))
                                DAmpz_rel = np.reshape(DAmpz_rel,np.size(DAmpz_rel))
                                DAmpx_rel_all = np.concatenate((DAmpx_rel_all,DAmpx_rel))
                                DAmpy_rel_all = np.concatenate((DAmpy_rel_all,DAmpy_rel))
                                DAmpz_rel_all = np.concatenate((DAmpz_rel_all,DAmpz_rel))

                                DAmpx = np.reshape(DAmpx,np.size(DAmpx))
                                DAmpy = np.reshape(DAmpy,np.size(DAmpy))
                                DAmpz = np.reshape(DAmpz,np.size(DAmpz))
                                DAmpx_all = np.concatenate((DAmpx_all,DAmpx))
                                DAmpy_all = np.concatenate((DAmpy_all,DAmpy))
                                DAmpz_all = np.concatenate((DAmpz_all,DAmpz))

                            else:
                                print('no antenna')
                                pass


    DAmpx_rel_all = [x for x in DAmpx_rel_all if not np.isnan(x)]
    DAmpy_rel_all = [x for x in DAmpy_rel_all if not np.isnan(x)]
    DAmpz_rel_all = [x for x in DAmpz_rel_all if not np.isnan(x)]
    DAmpx_all = [x for x in DAmpx_all if not np.isnan(x)]
    DAmpy_all = [x for x in DAmpy_all if not np.isnan(x)]
    DAmpz_all = [x for x in DAmpz_all if not np.isnan(x)]

    DAmpx_rel_min = np.amin(DAmpx_rel_all)
    DAmpx_rel_max = np.amax(DAmpx_rel_all)
    DAmpy_rel_min = np.amin(DAmpy_rel_all)
    DAmpy_rel_max = np.amax(DAmpy_rel_all)
    DAmpz_rel_min = np.amin(DAmpz_rel_all)
    DAmpz_rel_max = np.amax(DAmpz_rel_all)
    DAmpx_min = np.amin(DAmpx_all)
    DAmpx_max = np.amax(DAmpx_all)
    DAmpy_min = np.amin(DAmpy_all)
    DAmpy_max = np.amax(DAmpy_all)
    DAmpz_min = np.amin(DAmpz_all)
    DAmpz_max = np.amax(DAmpz_all)

    DAmpx_rel_bin = np.arange(DAmpx_rel_min,DAmpx_rel_max,(DAmpx_rel_max-DAmpx_rel_min)/1000.)
    DAmpy_rel_bin = np.arange(DAmpy_rel_min,DAmpy_rel_max,(DAmpy_rel_max-DAmpy_rel_min)/1000.)
    DAmpz_rel_bin = np.arange(DAmpz_rel_min,DAmpz_rel_max,(DAmpz_rel_max-DAmpz_rel_min)/1000.)
    DAmpx_bin = np.arange(DAmpx_min,DAmpx_max,(DAmpx_max-DAmpx_min)/1000.)
    DAmpy_bin = np.arange(DAmpy_min,DAmpy_max,(DAmpy_max-DAmpy_min)/1000.)
    DAmpz_bin = np.arange(DAmpz_min,DAmpz_max,(DAmpz_max-DAmpz_min)/1000.)

    '''
    #Gaussian
    mu_DAmpx_rel_all,std_DAmpx_rel_all = stats.norm.fit(DAmpx_rel_all)
    mu_DAmpy_rel_all,std_DAmpy_rel_all = stats.norm.fit(DAmpy_rel_all)
    mu_DAmpz_rel_all,std_DAmpz_rel_all = stats.norm.fit(DAmpz_rel_all)
    mu_DAmpx_all,std_DAmpx_all = stats.norm.fit(DAmpx_all)
    mu_DAmpy_all,std_DAmpy_all = stats.norm.fit(DAmpy_all)
    mu_DAmpz_all,std_DAmpz_all = stats.norm.fit(DAmpz_all)

    hist_fit_Dampx_rel = stats.norm.pdf(DAmpx_rel_bin,loc=mu_DAmpx_rel_all,scale=std_DAmpx_rel_all)
    hist_fit_Dampy_rel = stats.norm.pdf(DAmpy_rel_bin,loc=mu_DAmpy_rel_all,scale=std_DAmpy_rel_all)
    hist_fit_Dampz_rel = stats.norm.pdf(DAmpz_rel_bin,loc=mu_DAmpz_rel_all,scale=std_DAmpz_rel_all)
    hist_fit_Dampx = stats.norm.pdf(DAmpx_bin,loc=mu_DAmpx_all,scale=std_DAmpx_all)
    hist_fit_Dampy = stats.norm.pdf(DAmpy_bin,loc=mu_DAmpy_all,scale=std_DAmpy_all)
    hist_fit_Dampz = stats.norm.pdf(DAmpz_bin,loc=mu_DAmpz_all,scale=std_DAmpz_all)

    label_rel_x = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(mu_DAmpx_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpx_rel_all,2))
    label_rel_y = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(mu_DAmpy_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpy_rel_all,2))
    label_rel_z = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(mu_DAmpz_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpz_rel_all,2))

    labelx = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpx_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpx_all,2))
    labely = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpy_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpy_all,2))
    labelz = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpz_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpz_all,2))
    '''

    '''
    #Asymmetric Gaussian
    valskew_rel_x = stats.skew(DAmpx_rel_all)
    valskew_rel_y = stats.skew(DAmpx_rel_all)
    valskew_rel_z = stats.skew(DAmpx_rel_all)
    valskew_x = stats.skew(DAmpx_all)
    valskew_y = stats.skew(DAmpy_all)
    valskew_z = stats.skew(DAmpz_all)
    alpha_DAmpx_rel_all,loc_DAmpx_rel_all,scale_DAmpx_rel_all = stats.skewnorm.fit(DAmpx_rel_all)#,valskew_rel_x)
    alpha_DAmpy_rel_all,loc_DAmpy_rel_all,scale_DAmpy_rel_all = stats.skewnorm.fit(DAmpy_rel_all)#,valskew_rel_y)
    alpha_DAmpz_rel_all,loc_DAmpz_rel_all,scale_DAmpz_rel_all = stats.skewnorm.fit(DAmpz_rel_all)#,valskew_rel_z)
    alpha_DAmpx_all,loc_DAmpx_all,scale_DAmpx_all = stats.skewnorm.fit(DAmpx_all)#,valskew_x)
    alpha_DAmpy_all,loc_DAmpy_all,scale_DAmpy_all = stats.skewnorm.fit(DAmpy_all)#,valskew_y)
    alpha_DAmpz_all,loc_DAmpz_all,scale_DAmpz_all = stats.skewnorm.fit(DAmpz_all)#,valskew_z)

    hist_fit_Dampx_rel = stats.skewnorm.pdf(DAmpx_rel_bin,alpha_DAmpx_rel_all,loc=loc_DAmpx_rel_all,scale=scale_DAmpx_rel_all)
    hist_fit_Dampy_rel = stats.skewnorm.pdf(DAmpy_rel_bin,alpha_DAmpy_rel_all,loc=loc_DAmpy_rel_all,scale=scale_DAmpy_rel_all)
    hist_fit_Dampz_rel = stats.skewnorm.pdf(DAmpz_rel_bin,alpha_DAmpz_rel_all,loc=loc_DAmpz_rel_all,scale=scale_DAmpz_rel_all)
    hist_fit_Dampx = stats.skewnorm.pdf(DAmpx_bin,alpha_DAmpx_all,loc=loc_DAmpx_all,scale=scale_DAmpx_all)
    hist_fit_Dampy = stats.skewnorm.pdf(DAmpy_bin,alpha_DAmpy_all,loc=loc_DAmpy_all,scale=scale_DAmpy_all)
    hist_fit_Dampz = stats.skewnorm.pdf(DAmpz_bin,alpha_DAmpz_all,loc=loc_DAmpz_all,scale=scale_DAmpz_all)

    label_rel_x = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpx_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpx_rel_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpx_rel_all,2))
    label_rel_y = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpy_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpy_rel_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpy_rel_all,2))
    label_rel_z = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpz_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpz_rel_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpz_rel_all,2))

    labelx = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpx_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpx_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpx_all,2))
    labely = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpy_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpy_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpy_all,2))
    labelz = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpz_all,2))+'\n'+r' $\sigma$ = '+str(np.round(scale_DAmpz_all,2))+'\n'+r' $\alpha$ = '+str(np.round(alpha_DAmpz_all,2))
    '''

    '''
    #Gamma Law
    alpha_DAmpy_rel_all,loc_DAmpy_rel_all,scale_DAmpy_rel_all = stats.gamma.fit(DAmpy_rel_all)
    alpha_DAmpz_rel_all,loc_DAmpz_rel_all,scale_DAmpz_rel_all = stats.gamma.fit(DAmpz_rel_all)
    alpha_DAmpx_rel_all,loc_DAmpx_rel_all,scale_DAmpx_rel_all = stats.gamma.fit(DAmpx_rel_all)


    alpha_DAmpx_all,loc_DAmpx_all,scale_DAmpx_all = stats.gamma.fit(DAmpx_all)
    alpha_DAmpy_all,loc_DAmpy_all,scale_DAmpy_all = stats.gamma.fit(DAmpy_all)
    alpha_DAmpz_all,loc_DAmpz_all,scale_DAmpz_all = stats.gamma.fit(DAmpz_all)

    hist_fit_Dampx_rel = stats.gamma.pdf(DAmpx_rel_bin,alpha_DAmpx_rel_all,loc=loc_DAmpx_rel_all,scale=scale_DAmpx_rel_all)
    hist_fit_Dampy_rel = stats.gamma.pdf(DAmpy_rel_bin,alpha_DAmpy_rel_all,loc=loc_DAmpy_rel_all,scale=scale_DAmpy_rel_all)
    hist_fit_Dampz_rel = stats.gamma.pdf(DAmpz_rel_bin,alpha_DAmpz_rel_all,loc=loc_DAmpz_rel_all,scale=scale_DAmpz_rel_all)
    hist_fit_Dampx = stats.gamma.pdf(DAmpx_bin,alpha_DAmpx_all,loc=loc_DAmpx_all,scale=scale_DAmpx_all)
    hist_fit_Dampy = stats.gamma.pdf(DAmpy_bin,alpha_DAmpy_all,loc=loc_DAmpy_all,scale=scale_DAmpy_all)
    hist_fit_Dampz = stats.gamma.pdf(DAmpz_bin,alpha_DAmpz_all,loc=loc_DAmpz_all,scale=scale_DAmpz_all)

    label_rel_x = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpx_rel_all,2))+'\n'+r' $\sigma$ = '+str(scale_DAmpx_rel_all)
    label_rel_y = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpy_rel_all,2))+'\n'+r' $\sigma$ = '+str(scale_DAmpy_rel_all)
    label_rel_z = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(loc_DAmpz_rel_all,2))+'\n'+r' $\sigma$ = '+str(scale_DAmpz_rel_all)

    labelx = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\alpha = $'+str(np.round(loc_DAmpx_all,2))+'\n'+r' $\beta$ = '+str(np.round(scale_DAmpx_all,2))
    labely = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\alpha = $'+str(np.round(loc_DAmpy_all,2))+'\n'+r' $\beta$ = '+str(np.round(scale_DAmpy_all,2))
    labelz = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\alpha = $'+str(np.round(loc_DAmpz_all,2))+'\n'+r' $\beta$ = '+str(np.round(scale_DAmpz_all,2))
    '''

    '''
    #Gumbel right
    mu_DAmpx_rel_all,std_DAmpx_rel_all = stats.gumbel_r.fit(DAmpx_rel_all)
    mu_DAmpy_rel_all,std_DAmpy_rel_all = stats.gumbel_r.fit(DAmpy_rel_all)
    mu_DAmpz_rel_all,std_DAmpz_rel_all = stats.gumbel_r.fit(DAmpz_rel_all)
    mu_DAmpx_all,std_DAmpx_all = stats.gumbel_r.fit(DAmpx_all)
    mu_DAmpy_all,std_DAmpy_all = stats.gumbel_r.fit(DAmpy_all)
    mu_DAmpz_all,std_DAmpz_all = stats.gumbel_r.fit(DAmpz_all)

    hist_fit_Dampx_rel = stats.gumbel_r.pdf(DAmpx_rel_bin,loc=mu_DAmpx_rel_all,scale=std_DAmpx_rel_all)
    hist_fit_Dampy_rel = stats.gumbel_r.pdf(DAmpy_rel_bin,loc=mu_DAmpy_rel_all,scale=std_DAmpy_rel_all)
    hist_fit_Dampz_rel = stats.gumbel_r.pdf(DAmpz_rel_bin,loc=mu_DAmpz_rel_all,scale=std_DAmpz_rel_all)
    hist_fit_Dampx = stats.gumbel_r.pdf(DAmpx_bin,loc=mu_DAmpx_all,scale=std_DAmpx_all)
    hist_fit_Dampy = stats.gumbel_r.pdf(DAmpy_bin,loc=mu_DAmpy_all,scale=std_DAmpy_all)
    hist_fit_Dampz = stats.gumbel_r.pdf(DAmpz_bin,loc=mu_DAmpz_all,scale=std_DAmpz_all)

    label_rel_x = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpx_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpx_rel_all,2))
    label_rel_y = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpy_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpy_rel_all,2))
    label_rel_z = r'$\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpz_rel_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpz_rel_all,2))

    labelx = r'Fit $\rm \Delta E_{x}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpx_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpx_all,2))
    labely = r'Fit $\rm \Delta E_{y}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpy_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpy_all,2))
    labelz = r'Fit $\rm \Delta E_{z}$'+'\n'+r'$\mu = $'+str(np.round(mu_DAmpz_all,2))+'\n'+r' $\sigma$ = '+str(np.round(std_DAmpz_all,2))
    '''

    nbin='auto' #30
    fig1 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    print(np.shape(DAmpx_rel_all))
    hist_DAmpx_rel_all,bin_DAmpx_rel_all,tmp = pl.hist(DAmpx_rel_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_rel_bin,hist_fit_Dampx_rel,'--r',label=label_rel_x,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{x}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_rel_all,bin_DAmpy_rel_all,tmp = pl.hist(DAmpy_rel_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_rel_bin,hist_fit_Dampy_rel,'--r',label=label_rel_y,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{y}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_rel_all,bin_DAmpz_rel_all,tmp = pl.hist(DAmpz_rel_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_rel_bin,hist_fit_Dampz_rel,'--r',label=label_rel_z,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{z}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_rel_all,DAmpy_rel_all,DAmpz_rel_all],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Relative difference on E')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/RelDeltaE_high.png'
    pl.savefig(figname,dpi=350)

    fig2 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    hist_DAmpx_all,bin_DAmpx_all,tmp = pl.hist(DAmpx_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_bin,hist_fit_Dampx,'--r',label=labelx,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{x}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_all,bin_DAmpy_all,tmp = pl.hist(DAmpy_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_bin,hist_fit_Dampy,'--r',label=labely,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{y}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_all,bin_DAmpz_all,tmp = pl.hist(DAmpz_all,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_bin,hist_fit_Dampz,'--r',label=labelz,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{z}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_all,DAmpy_all,DAmpz_all],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Absolute difference on E')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/AbsDeltaE.png'
    pl.savefig(figname,dpi=350)

    print('\n mu Rel DEx',np.round(np.mean(DAmpx_rel_all),2),'mu Rel DEy',np.round(np.mean(DAmpy_rel_all),2),'mu Rel DEz',np.round(np.mean(DAmpz_rel_all),2))
    print(' std Rel DEx',np.round(np.std(DAmpx_rel_all),2),'std Rel DEy',np.round(np.std(DAmpy_rel_all),2),'std Rel DEz',np.round(np.std(DAmpz_rel_all),2))
    print(' med Rel DEx',np.round(np.median(DAmpx_rel_all),2),'med Rel DEy',np.round(np.median(DAmpy_rel_all),2),'med Rel DEz',np.round(np.median(DAmpz_rel_all),2))
    print(' pos max Rel DEx',np.round(bin_DAmpx_rel_all[np.argmax(hist_DAmpx_rel_all)],2),'pos max Rel DEy',np.round(bin_DAmpy_rel_all[np.argmax(hist_DAmpy_rel_all)],2),'pos max Rel DEz',np.round(bin_DAmpz_rel_all[np.argmax(hist_DAmpz_rel_all)],2))
    print('\n mu Abso DEx',np.round(np.mean(DAmpx_all),2),'mu Abso DEy',np.round(np.mean(DAmpy_all),2),'mu Abso DEz',np.round(np.mean(DAmpz_all),2))
    print(' std Abso DEx',np.round(np.std(DAmpx_all),2),'std Abso DEy',np.round(np.std(DAmpy_all),2),'std Abso DEz',np.round(np.std(DAmpz_all),2))
    print(' med Abso DEx',np.round(np.median(DAmpx_all),2),'med Abso DEy',np.round(np.median(DAmpy_all),2),'med Abso DEz',np.round(np.median(DAmpz_all),2))
    print(' pos max Abso DEx',np.round(bin_DAmpx_all[np.argmax(hist_DAmpx_all)],2),'pos max Abso DEy',np.round(bin_DAmpy_all[np.argmax(hist_DAmpy_all)],2),'pos max Abso DEz',np.round(bin_DAmpz_all[np.argmax(hist_DAmpz_all)],2))
    print(' ')

    if DISPLAY:
        pl.show()
    pl.close(fig1)
    pl.close(fig2)

    return

############################################################################################################
def histo_DEpp_threshold(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires):

    nbin='auto' #30
    threshold = 0.4
    DAmpx_rel_all_sup = []
    DAmpy_rel_all_sup = []
    DAmpz_rel_all_sup = []
    DAmpx_all_sup = []
    DAmpy_all_sup = []
    DAmpz_all_sup = []
    DAmpx_rel_all_inf = []
    DAmpy_rel_all_inf = []
    DAmpz_rel_all_inf = []
    DAmpx_all_inf = []
    DAmpy_all_inf = []
    DAmpz_all_inf = []
    for iserie in Nserie:
        for ie in range(len(Etarget)):
            ratioE = float(Etarget[ie])/1E+19

            for iDd in Dd:
                for ialfa in range(np.size(alfa)):
                    if alfa[ialfa]>0 or (alfa[ialfa]==0 and iDd==20000):
                        folder_zhaires=ROOT_zhaires+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+'/'
                        folder_RM=ROOT_RM+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'/'
                        showers = glob.glob(folder_zhaires+'[0-9][0-9]*[0-9]/')
                        #showers=showers[0:11]

                        #Start loop on shower
                        for shower in showers:
                            showerID = shower.replace('/','.').split('.')[-2]
                            trace_files=glob.glob(shower+'/split/*.trace')
                            Nant=len(trace_files)
                            if Nant>0:
                                Ampx_RM,Ampy_RM,Ampz_RM,x0_RM,y0_RM,z0_RM = Efield_RM(folder_RM+showerID+'/')
                                Ampx_zhaires,Ampy_zhaires,Ampz_zhaires,x0_zhaires,y0_zhaires,z0_zhaires = Efield_zhaires(shower)

                                ind = np.where(Ampy_RM==0)

                                if np.shape(ind)[1]!=np.size(Ampy_RM):
                                    Ampx_RM = np.delete(Ampx_RM,ind)
                                    Ampy_RM = np.delete(Ampy_RM,ind)
                                    Ampz_RM = np.delete(Ampz_RM,ind)
                                    Ampx_zhaires = np.delete(Ampx_zhaires,ind)
                                    Ampy_zhaires = np.delete(Ampy_zhaires,ind)
                                    Ampz_zhaires = np.delete(Ampz_zhaires,ind)

                                    # Split the arrays in two: high and low voltage
                                    indx_sup = np.where(Ampx_zhaires>threshold*np.amax(Ampx_zhaires))
                                    indy_sup = np.where(Ampy_zhaires>threshold*np.amax(Ampy_zhaires))
                                    indz_sup = np.where(Ampz_zhaires>threshold*np.amax(Ampz_zhaires))
                                    indx_inf = np.where(Ampx_zhaires<=threshold*np.amax(Ampx_zhaires))
                                    indy_inf = np.where(Ampy_zhaires<=threshold*np.amax(Ampy_zhaires))
                                    indz_inf = np.where(Ampz_zhaires<=threshold*np.amax(Ampz_zhaires))

                                    Ampx_RM_sup = Ampx_RM[indx_sup]
                                    Ampy_RM_sup = Ampy_RM[indy_sup]
                                    Ampz_RM_sup = Ampz_RM[indz_sup]
                                    Ampx_zhaires_sup = Ampx_zhaires[indx_sup]
                                    Ampy_zhaires_sup = Ampy_zhaires[indy_sup]
                                    Ampz_zhaires_sup = Ampz_zhaires[indz_sup]
                                    Ampx_RM_inf = Ampx_RM[indx_inf]
                                    Ampy_RM_inf = Ampy_RM[indy_inf]
                                    Ampz_RM_inf = Ampz_RM[indz_inf]
                                    Ampx_zhaires_inf = Ampx_zhaires[indx_inf]
                                    Ampy_zhaires_inf = Ampy_zhaires[indy_inf]
                                    Ampz_zhaires_inf = Ampz_zhaires[indz_inf]

                                    ## Computation part
                                    DAmpx_sup = Ampx_RM_sup-Ampx_zhaires_sup
                                    DAmpy_sup = Ampy_RM_sup-Ampy_zhaires_sup
                                    DAmpz_sup = Ampz_RM_sup-Ampz_zhaires_sup
                                    DAmpx_rel_sup = (Ampx_RM_sup-Ampx_zhaires_sup)/Ampx_zhaires_sup
                                    DAmpy_rel_sup = (Ampy_RM_sup-Ampy_zhaires_sup)/Ampy_zhaires_sup
                                    DAmpz_rel_sup = (Ampz_RM_sup-Ampz_zhaires_sup)/Ampz_zhaires_sup
                                    DAmpx_inf = Ampx_RM_inf-Ampx_zhaires_inf
                                    DAmpy_inf = Ampy_RM_inf-Ampy_zhaires_inf
                                    DAmpz_inf = Ampz_RM_inf-Ampz_zhaires_inf
                                    DAmpx_rel_inf = (Ampx_RM_inf-Ampx_zhaires_inf)/Ampx_zhaires_inf
                                    DAmpy_rel_inf = (Ampy_RM_inf-Ampy_zhaires_inf)/Ampy_zhaires_inf
                                    DAmpz_rel_inf = (Ampz_RM_inf-Ampz_zhaires_inf)/Ampz_zhaires_inf

                                    DAmpx_rel_sup = np.reshape(DAmpx_rel_sup,np.size(DAmpx_rel_sup))
                                    DAmpy_rel_sup = np.reshape(DAmpy_rel_sup,np.size(DAmpy_rel_sup))
                                    DAmpz_rel_sup = np.reshape(DAmpz_rel_sup,np.size(DAmpz_rel_sup))
                                    DAmpx_rel_all_sup = np.concatenate((DAmpx_rel_all_sup,DAmpx_rel_sup))
                                    DAmpy_rel_all_sup = np.concatenate((DAmpy_rel_all_sup,DAmpy_rel_sup))
                                    DAmpz_rel_all_sup = np.concatenate((DAmpz_rel_all_sup,DAmpz_rel_sup))
                                    DAmpx_rel_inf = np.reshape(DAmpx_rel_inf,np.size(DAmpx_rel_inf))
                                    DAmpy_rel_inf = np.reshape(DAmpy_rel_inf,np.size(DAmpy_rel_inf))
                                    DAmpz_rel_inf = np.reshape(DAmpz_rel_inf,np.size(DAmpz_rel_inf))
                                    DAmpx_rel_all_inf = np.concatenate((DAmpx_rel_all_inf,DAmpx_rel_inf))
                                    DAmpy_rel_all_inf = np.concatenate((DAmpy_rel_all_inf,DAmpy_rel_inf))
                                    DAmpz_rel_all_inf = np.concatenate((DAmpz_rel_all_inf,DAmpz_rel_inf))

                                    DAmpx_sup = np.reshape(DAmpx_sup,np.size(DAmpx_sup))
                                    DAmpy_sup = np.reshape(DAmpy_sup,np.size(DAmpy_sup))
                                    DAmpz_sup = np.reshape(DAmpz_sup,np.size(DAmpz_sup))
                                    DAmpx_all_sup = np.concatenate((DAmpx_all_sup,DAmpx_sup))
                                    DAmpy_all_sup = np.concatenate((DAmpy_all_sup,DAmpy_sup))
                                    DAmpz_all_sup = np.concatenate((DAmpz_all_sup,DAmpz_sup))
                                    DAmpx_inf = np.reshape(DAmpx_inf,np.size(DAmpx_inf))
                                    DAmpy_inf = np.reshape(DAmpy_inf,np.size(DAmpy_inf))
                                    DAmpz_inf = np.reshape(DAmpz_inf,np.size(DAmpz_inf))
                                    DAmpx_all_inf = np.concatenate((DAmpx_all_inf,DAmpx_inf))
                                    DAmpy_all_inf = np.concatenate((DAmpy_all_inf,DAmpy_inf))
                                    DAmpz_all_inf = np.concatenate((DAmpz_all_inf,DAmpz_inf))

                            else:
                                print('no antenna')
                                pass

    DAmpx_rel_all_sup = [x for x in DAmpx_rel_all_sup if not np.isnan(x)]
    DAmpy_rel_all_sup = [x for x in DAmpy_rel_all_sup if not np.isnan(x)]
    DAmpz_rel_all_sup = [x for x in DAmpz_rel_all_sup if not np.isnan(x)]
    DAmpx_all_sup = [x for x in DAmpx_all_sup if not np.isnan(x)]
    DAmpy_all_sup = [x for x in DAmpy_all_sup if not np.isnan(x)]
    DAmpz_all_sup = [x for x in DAmpz_all_sup if not np.isnan(x)]
    DAmpx_rel_all_inf = [x for x in DAmpx_rel_all_inf if not np.isnan(x)]
    DAmpy_rel_all_inf = [x for x in DAmpy_rel_all_inf if not np.isnan(x)]
    DAmpz_rel_all_inf = [x for x in DAmpz_rel_all_inf if not np.isnan(x)]
    DAmpx_all_inf = [x for x in DAmpx_all_inf if not np.isnan(x)]
    DAmpy_all_inf = [x for x in DAmpy_all_inf if not np.isnan(x)]
    DAmpz_all_inf = [x for x in DAmpz_all_inf if not np.isnan(x)]

    DAmpx_rel_min_sup = np.amin(DAmpx_rel_all_sup)
    DAmpx_rel_max_sup = np.amax(DAmpx_rel_all_sup)
    DAmpy_rel_min_sup = np.amin(DAmpy_rel_all_sup)
    DAmpy_rel_max_sup = np.amax(DAmpy_rel_all_sup)
    DAmpz_rel_min_sup = np.amin(DAmpz_rel_all_sup)
    DAmpz_rel_max_sup = np.amax(DAmpz_rel_all_sup)
    DAmpx_min_sup = np.amin(DAmpx_all_sup)
    DAmpx_max_sup = np.amax(DAmpx_all_sup)
    DAmpy_min_sup = np.amin(DAmpy_all_sup)
    DAmpy_max_sup = np.amax(DAmpy_all_sup)
    DAmpz_min_sup = np.amin(DAmpz_all_sup)
    DAmpz_max_sup = np.amax(DAmpz_all_sup)
    DAmpx_rel_min_inf = np.amin(DAmpx_rel_all_inf)
    DAmpx_rel_max_inf = np.amax(DAmpx_rel_all_inf)
    DAmpy_rel_min_inf = np.amin(DAmpy_rel_all_inf)
    DAmpy_rel_max_inf = np.amax(DAmpy_rel_all_inf)
    DAmpz_rel_min_inf = np.amin(DAmpz_rel_all_inf)
    DAmpz_rel_max_inf = np.amax(DAmpz_rel_all_inf)
    DAmpx_min_inf = np.amin(DAmpx_all_inf)
    DAmpx_max_inf = np.amax(DAmpx_all_inf)
    DAmpy_min_inf = np.amin(DAmpy_all_inf)
    DAmpy_max_inf = np.amax(DAmpy_all_inf)
    DAmpz_min_inf = np.amin(DAmpz_all_inf)
    DAmpz_max_inf = np.amax(DAmpz_all_inf)

    DAmpx_rel_bin_sup = np.arange(DAmpx_rel_min_sup,DAmpx_rel_max_sup,(DAmpx_rel_max_sup-DAmpx_rel_min_sup)/1000.)
    DAmpy_rel_bin_sup = np.arange(DAmpy_rel_min_sup,DAmpy_rel_max_sup,(DAmpy_rel_max_sup-DAmpy_rel_min_sup)/1000.)
    DAmpz_rel_bin_sup = np.arange(DAmpz_rel_min_sup,DAmpz_rel_max_sup,(DAmpz_rel_max_sup-DAmpz_rel_min_sup)/1000.)
    DAmpx_bin_sup = np.arange(DAmpx_min_sup,DAmpx_max_sup,(DAmpx_max_sup-DAmpx_min_sup)/1000.)
    DAmpy_bin_sup = np.arange(DAmpy_min_sup,DAmpy_max_sup,(DAmpy_max_sup-DAmpy_min_sup)/1000.)
    DAmpz_bin_sup = np.arange(DAmpz_min_sup,DAmpz_max_sup,(DAmpz_max_sup-DAmpz_min_sup)/1000.)
    DAmpx_rel_bin_inf = np.arange(DAmpx_rel_min_inf,DAmpx_rel_max_inf,(DAmpx_rel_max_inf-DAmpx_rel_min_inf)/1000.)
    DAmpy_rel_bin_inf = np.arange(DAmpy_rel_min_inf,DAmpy_rel_max_inf,(DAmpy_rel_max_inf-DAmpy_rel_min_inf)/1000.)
    DAmpz_rel_bin_inf = np.arange(DAmpz_rel_min_inf,DAmpz_rel_max_inf,(DAmpz_rel_max_inf-DAmpz_rel_min_inf)/1000.)
    DAmpx_bin_inf = np.arange(DAmpx_min_inf,DAmpx_max_inf,(DAmpx_max_inf-DAmpx_min_inf)/1000.)
    DAmpy_bin_inf = np.arange(DAmpy_min_inf,DAmpy_max_inf,(DAmpy_max_inf-DAmpy_min_inf)/1000.)
    DAmpz_bin_inf = np.arange(DAmpz_min_inf,DAmpz_max_inf,(DAmpz_max_inf-DAmpz_min_inf)/1000.)

    #>threshold
    fig1 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    hist_DAmpx_rel_all_sup,bin_DAmpx_rel_all_sup,tmp = pl.hist(DAmpx_rel_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_rel_bin_sup,hist_fit_Dampx_rel_sup,'--r',label=label_rel_x,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{x}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_rel_all_sup,bin_DAmpy_rel_all_sup,tmp = pl.hist(DAmpy_rel_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_rel_bin_sup,hist_fit_Dampy_rel_sup,'--r',label=label_rel_y,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{y}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_rel_all_sup,bin_DAmpz_rel_all_sup,tmp = pl.hist(DAmpz_rel_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_rel_bin_sup,hist_fit_Dampz_rel_sup,'--r',label=label_rel_z,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{z}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_rel_all_sup,DAmpy_rel_all_sup,DAmpz_rel_all_sup],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Relative difference on E$>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/RelDeltaE_high.png'
    pl.savefig(figname,dpi=350)

    fig2 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    hist_DAmpx_all_sup,bin_DAmpx_all_sup,tmp = pl.hist(DAmpx_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_bin_sup,hist_fit_Dampx_sup,'--r',label=labelx,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{x}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_all_sup,bin_DAmpy_all_sup,tmp = pl.hist(DAmpy_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_bin_sup,hist_fit_Dampy_sup,'--r',label=labely,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{y}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_all_sup,bin_DAmpz_all_sup,tmp = pl.hist(DAmpz_all_sup,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_bin_sup,hist_fit_Dampz_sup,'--r',label=labelz,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{z}>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_all_sup,DAmpy_rel_all_sup,DAmpz_rel_all_sup],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Absolute difference on E$>$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/AbsDeltaE_high.png'
    pl.savefig(figname,dpi=350)

    print('\n mu Rel DEx',np.round(np.mean(DAmpx_rel_all_sup),2),'mu Rel DEy',np.round(np.mean(DAmpy_rel_all_sup),2),'mu Rel DEz',np.round(np.mean(DAmpz_rel_all_sup),2))
    print(' std Rel DEx',np.round(np.std(DAmpx_rel_all_sup),2),'std Rel DEy',np.round(np.std(DAmpy_rel_all_sup),2),'std Rel DEz',np.round(np.std(DAmpz_rel_all_sup),2))
    print(' med Rel DEx',np.round(np.median(DAmpx_rel_all_sup),2),'med Rel DEy',np.round(np.median(DAmpy_rel_all_sup),2),'med Rel DEz',np.round(np.median(DAmpz_rel_all_sup),2))
    print(' pos max Rel DEx',np.round(bin_DAmpx_rel_all_sup[np.argmax(hist_DAmpx_rel_all_sup)],2),'pos max Rel DEy',np.round(bin_DAmpy_rel_all_sup[np.argmax(hist_DAmpy_rel_all_sup)],2),'pos max Rel DEz',np.round(bin_DAmpz_rel_all_sup[np.argmax(hist_DAmpz_rel_all_sup)],2))
    print('\n mu Abso DEx',np.round(np.mean(DAmpx_all_sup),2),'mu Abso DEy',np.round(np.mean(DAmpy_all_sup),2),'mu Abso DEz',np.round(np.mean(DAmpz_all_sup),2))
    print(' std Abso DEx',np.round(np.std(DAmpx_all_sup),2),'std Abso DEy',np.round(np.std(DAmpy_all_sup),2),'std Abso DEz',np.round(np.std(DAmpz_all_sup),2))
    print(' med Abso DEx',np.round(np.median(DAmpx_all_sup),2),'med Abso DEy',np.round(np.median(DAmpy_all_sup),2),'med Abso DEz',np.round(np.median(DAmpz_all_sup),2))
    print(' pos max Abso DEx',np.round(bin_DAmpx_all_sup[np.argmax(hist_DAmpx_all_sup)],2),'pos max Abso DEy',np.round(bin_DAmpy_all_sup[np.argmax(hist_DAmpy_all_sup)],2),'pos max Abso DEz',np.round(bin_DAmpz_all_sup[np.argmax(hist_DAmpz_all_sup)],2))
    print(' ')

    ############
    #<threshold
    fig3 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    hist_DAmpx_rel_all_inf,bin_DAmpx_rel_all_inf,tmp = pl.hist(DAmpx_rel_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_rel_bin,hist_fit_Dampx_rel,'--r',label=label_rel_x,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{x}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_rel_all_inf,bin_DAmpy_rel_all_inf,tmp = pl.hist(DAmpy_rel_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_rel_bin,hist_fit_Dampy_rel,'--r',label=label_rel_y,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{y}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_rel_all_inf,bin_DAmpz_rel_all_inf,tmp = pl.hist(DAmpz_rel_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_rel_bin,hist_fit_Dampz_rel,'--r',label=label_rel_z,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{z}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_rel_all_inf,DAmpy_rel_all_inf,DAmpz_rel_all_inf],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Relative difference on E$<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/RelDeltaE_low.png'
    pl.savefig(figname,dpi=350)

    fig4 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    hist_DAmpx_all_inf,bin_DAmpx_all_inf,tmp = pl.hist(DAmpx_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpx_bin,hist_fit_Dampx,'--r',label=labelx,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{x}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    hist_DAmpy_all_inf,bin_DAmpy_all_inf,tmp = pl.hist(DAmpy_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpy_bin,hist_fit_Dampy,'--r',label=labely,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{y}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    hist_DAmpz_all_inf,bin_DAmpz_all_inf,tmp = pl.hist(DAmpz_all_inf,bins=nbin,normed=1)
    try:
        pl.plot(DAmpz_bin,hist_fit_Dampz,'--r',label=labelz,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{z}<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.subplot(224)
    pl.hist([DAmpx_all_inf,DAmpy_rel_all_inf,DAmpz_rel_all_inf],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    pl.xlabel('Absolute difference on E$<$'+str(threshold)+'E$_{max}$')
    pl.ylabel('Rate')
    pl.legend()
    figname = figdir+'/AbsDeltaE_low.png'
    pl.savefig(figname,dpi=350)

    print('\n mu Rel DEx',np.round(np.mean(DAmpx_rel_all_inf),2),'mu Rel DEy',np.round(np.mean(DAmpy_rel_all_inf),2),'mu Rel DEz',np.round(np.mean(DAmpz_rel_all_inf),2))
    print(' std Rel DEx',np.round(np.std(DAmpx_rel_all_inf),2),'std Rel DEy',np.round(np.std(DAmpy_rel_all_inf),2),'std Rel DEz',np.round(np.std(DAmpz_rel_all_inf),2))
    print(' med Rel DEx',np.round(np.median(DAmpx_rel_all_inf),2),'med Rel DEy',np.round(np.median(DAmpy_rel_all_inf),2),'med Rel DEz',np.round(np.median(DAmpz_rel_all_inf),2))
    print(' pos max Rel DEx',np.round(bin_DAmpx_rel_all_inf[np.argmax(hist_DAmpx_rel_all_inf)],2),'pos max Rel DEy',np.round(bin_DAmpy_rel_all_inf[np.argmax(hist_DAmpy_rel_all_inf)],2),'pos max Rel DEz',np.round(bin_DAmpz_rel_all_inf[np.argmax(hist_DAmpz_rel_all_inf)],2))
    print('\n mu Abso DEx',np.round(np.mean(DAmpx_all_inf),2),'mu Abso DEy',np.round(np.mean(DAmpy_all_inf),2),'mu Abso DEz',np.round(np.mean(DAmpz_all_inf),2))
    print(' std Abso DEx',np.round(np.std(DAmpx_all_inf),2),'std Abso DEy',np.round(np.std(DAmpy_all_inf),2),'std Abso DEz',np.round(np.std(DAmpz_all_inf),2))
    print(' med Abso DEx',np.round(np.median(DAmpx_all_inf),2),'med Abso DEy',np.round(np.median(DAmpy_all_inf),2),'med Abso DEz',np.round(np.median(DAmpz_all_inf),2))
    print(' pos max Abso DEx',np.round(bin_DAmpx_all_inf[np.argmax(hist_DAmpx_all_inf)],2),'pos max Abso DEy',np.round(bin_DAmpy_all_inf[np.argmax(hist_DAmpy_all_inf)],2),'pos max Abso DEz',np.round(bin_DAmpz_all_inf[np.argmax(hist_DAmpz_all_inf)],2))
    print(' ')

    if DISPLAY:
        pl.show()
    pl.close(fig1)
    pl.close(fig2)
    pl.close(fig3)
    pl.close(fig4)

    return

############################################################################################################
def histo_Ntrig(Dd,alfa,hz,Etarget,sep,ROOT_RM,ROOT_zhaires):
    task = 'EE'+Etarget[0]+'_D'+str(Dd[0])+'_alpha'+str(alfa[0])+'_height'+str(hz)+'_sep'+str(sep[0])
    det_file_zhaires=ROOT_zhaires+'detection/detection_count_'+task+'_serie1.txt'
    det_file_RM=ROOT_RM+'detection/detection_count_'+task+'.txt'
    res_zhaires = np.loadtxt(det_file_zhaires)
    res_RM = np.loadtxt(det_file_RM)

    Nant_ew_cons_zhaires = res_zhaires[:,3]
    Nant_ns_cons_zhaires = res_zhaires[:,4]
    Nant_up_cons_zhaires = res_zhaires[:,5]
    Nant_tot_cons_zhaires = res_zhaires[:,6]
    Nant_ew_agg_zhaires = res_zhaires[:,7]
    Nant_ns_agg_zhaires = res_zhaires[:,8]
    Nant_up_agg_zhaires = res_zhaires[:,9]
    Nant_tot_agg_zhaires = res_zhaires[:,10]
    Vp2p_ew_zhaires = res_zhaires[:,11]
    Vp2p_ns_zhaires = res_zhaires[:,12]
    Vp2p_up_zhaires = res_zhaires[:,13]
    Vp2p_tot_zhaires = res_zhaires[:,14]

    Nant_ew_cons_RM = res_RM[:,3]
    Nant_ns_cons_RM = res_RM[:,4]
    Nant_up_cons_RM = res_RM[:,5]
    Nant_tot_cons_RM = res_RM[:,6]
    Nant_ew_agg_RM = res_RM[:,7]
    Nant_ns_agg_RM = res_RM[:,8]
    Nant_up_agg_RM = res_RM[:,9]
    Nant_tot_agg_RM = res_RM[:,10]

    Vp2p_ew_RM = res_RM[:,11]
    Vp2p_ns_RM = res_RM[:,12]
    Vp2p_up_RM = res_RM[:,13]
    Vp2p_tot_RM = res_RM[:,14]

    nbin = 20
    fig1 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(241)
    pl.hist([Nant_ew_cons_zhaires,Nant_ew_cons_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{EW}$ Conservative')
    pl.subplot(242)
    pl.hist([Nant_ew_agg_zhaires,Nant_ew_agg_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{EW}$ Aggressive')
    pl.subplot(243)
    pl.hist([Nant_ns_cons_zhaires,Nant_ns_cons_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{NS}$ Conservative')
    pl.subplot(244)
    pl.hist([Nant_ns_agg_zhaires,Nant_ns_agg_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('WeiRateght')
    pl.legend()
    pl.title('V$_{NS}$ Aggressive')
    pl.subplot(245)
    pl.hist([Nant_up_cons_zhaires,Nant_up_cons_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{up}$ Conservative')
    pl.subplot(246)
    pl.hist([Nant_up_agg_zhaires,Nant_up_agg_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{up}$ Aggressive')
    pl.subplot(247)
    pl.hist([Nant_tot_cons_zhaires,Nant_tot_cons_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{tot}$ Conservative')
    pl.subplot(248)
    pl.hist([Nant_tot_agg_zhaires,Nant_tot_agg_RM],bins=nbin,normed=1,label=['ZHAireS','RadioMorphing'])
    pl.xlabel('Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('V$_{Tot}$ Aggressive')
    fig1.tight_layout()
    figname = figdir+'/Nant_Triggered.png'
    pl.savefig(figname,dpi=350)

    # Difference number of triggered antennas
    Dtrig_ew_cons = Nant_ew_cons_zhaires-Nant_ew_cons_RM
    Dtrig_ew_agg = Nant_ew_agg_zhaires-Nant_ew_agg_RM
    Dtrig_ns_cons = Nant_ns_cons_zhaires-Nant_ns_cons_RM
    Dtrig_ns_agg = Nant_ns_agg_zhaires-Nant_ns_agg_RM
    Dtrig_up_cons = Nant_up_cons_zhaires-Nant_up_cons_RM
    Dtrig_up_agg = Nant_up_agg_zhaires-Nant_up_agg_RM
    Dtrig_tot_cons = Nant_tot_cons_zhaires-Nant_tot_cons_RM
    Dtrig_tot_agg = Nant_tot_agg_zhaires-Nant_tot_agg_RM

    fig2 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(211)
    pl.hist([Dtrig_ew_cons,Dtrig_ns_cons,Dtrig_up_cons,Dtrig_tot_cons],bins=nbin,normed=1,label=['EW','NS','up','total'])
    pl.xlabel('$\Delta$ Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Conservative case')
    pl.subplot(212)
    pl.hist([Dtrig_ew_agg,Dtrig_ns_agg,Dtrig_up_agg,Dtrig_tot_agg],bins=nbin,normed=1,label=['EW','NS','up','total'])
    pl.xlabel('$\Delta$ Nant triggered')
    pl.ylabel('Rate')
    pl.legend()
    pl.title('Aggressive case')
    fig2.tight_layout()
    figname = figdir+'/DeltaNant.png'
    pl.savefig(figname,dpi=350)

    # Number of detected showers
    nantmin = 8
    ok_agg_RM = np.logical_or(Nant_ew_agg_RM>=nantmin,Nant_ns_agg_RM>=nantmin,Nant_up_agg_RM>=nantmin)
    ok_cons_RM = np.logical_or(Nant_ew_cons_RM>=nantmin,Nant_ns_cons_RM>=nantmin,Nant_up_cons_RM>=nantmin)
    ok_agg_zhaires = np.logical_or(Nant_ew_agg_zhaires>=nantmin,Nant_ns_agg_zhaires>=nantmin,Nant_up_agg_zhaires>=nantmin)
    ok_cons_zhaires = np.logical_or(Nant_ew_cons_zhaires>=nantmin,Nant_ns_cons_zhaires>=nantmin,Nant_up_cons_zhaires>=nantmin)

    not_ok_agg_RM = np.logical_not(ok_agg_RM)
    not_ok_agg_zhaires = np.logical_not(ok_agg_zhaires)

    print('Shower trigger for ZHAireS and RM:',np.sum((ok_agg_RM & ok_agg_zhaires)),', for RM but not ZHAireS:',np.sum((ok_agg_RM & not_ok_agg_zhaires)),', for ZHAireS but not RM:',np.sum((not_ok_agg_RM & ok_agg_zhaires)),', for none of them:',np.sum((not_ok_agg_RM & not_ok_agg_zhaires)))

    Nsh_trig_agg_RM = float(np.sum(ok_agg_RM))
    Nsh_trig_cons_RM = float(np.sum(ok_cons_RM))
    Nsh_trig_agg_zhaires = float(np.sum(ok_agg_zhaires))
    Nsh_trig_cons_zhaires = float(np.sum(ok_cons_zhaires))

    print('N triggered showers (Aggressive case) with RM:',Nsh_trig_agg_RM,' and with ZHAireS:',Nsh_trig_agg_zhaires)
    print('Absolute Delta triggered showers (Aggressive case):',Nsh_trig_agg_zhaires-Nsh_trig_agg_RM)
    print('Relative Delta triggered showers (Aggressive case):',100*(Nsh_trig_agg_zhaires-Nsh_trig_agg_RM)/Nsh_trig_agg_zhaires)
    print('N triggered showers (Conservative case) with RM:',Nsh_trig_cons_RM,' and with ZHAireS:',Nsh_trig_cons_zhaires)
    print('Absolute Delta triggered showers (Conservative case):',Nsh_trig_cons_zhaires-Nsh_trig_cons_RM)
    print('Relative Delta triggered showers (Conservative case):',100*(Nsh_trig_cons_zhaires-Nsh_trig_cons_RM)/Nsh_trig_cons_zhaires)
    print('Total number of showers:',np.size(ok_agg_zhaires))

    if DISPLAY:
        pl.show()
    pl.close(fig1)
    pl.close(fig2)

    return

############################################################################################################
def histo_DEpp_EnAzZen(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires):

    zen_bin,az_bin,eny_bin,injh_bin = build_binning()
    DAmpx_rel_all = [[]]*(len(zen_bin)+1)
    DAmpy_rel_all = [[]]*(len(zen_bin)+1)
    DAmpz_rel_all = [[]]*(len(zen_bin)+1)
    DAmpx_all = [[]]*(len(zen_bin)+1)
    DAmpy_all = [[]]*(len(zen_bin)+1)
    DAmpz_all = [[]]*(len(zen_bin)+1)
    for iserie in Nserie:
        for ie in range(len(Etarget)):
            ratioE = float(Etarget[ie])/1E+19

            for iDd in Dd:
                for ialfa in range(np.size(alfa)):
                    if alfa[ialfa]>0 or (alfa[ialfa]==0 and iDd==20000):
                        folder_zhaires=ROOT_zhaires+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+'/'
                        folder_RM=ROOT_RM+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'/'
                        showers = glob.glob(folder_zhaires+'[0-9][0-9]*[0-9]/')
                        showers=showers[0:11]

                        showerID_sh,zen_sh,az_sh,eny_sh,injh_sh = read_showers(showers)
                        zen_ind,az_ind,eny_ind,injh_ind = find_ind(showerID_sh,zen_sh,az_sh,eny_sh,injh_sh, zen_bin,az_bin,eny_bin,injh_bin)

                        #Start loop on shower
                        for shower in showers:
                            showerID = shower.replace('/','.').split('.')[-2]
                            trace_files=glob.glob(shower+'/split/*.trace')
                            Nant=len(trace_files)
                            if Nant>0:
                                Ampx_RM,Ampy_RM,Ampz_RM,x0_RM,y0_RM,z0_RM = Efield_RM(folder_RM+showerID+'/')
                                Ampx_zhaires,Ampy_zhaires,Ampz_zhaires,x0_zhaires,y0_zhaires,z0_zhaires = Efield_zhaires(shower)

                                ind = np.where(Ampy_RM==0)
                                Ampx_RM = np.delete(Ampx_RM,ind)
                                Ampy_RM = np.delete(Ampy_RM,ind)
                                Ampz_RM = np.delete(Ampz_RM,ind)
                                Ampx_zhaires = np.delete(Ampx_zhaires,ind)
                                Ampy_zhaires = np.delete(Ampy_zhaires,ind)
                                Ampz_zhaires = np.delete(Ampz_zhaires,ind)

                                ## Computation part
                                DAmpx = Ampx_RM-Ampx_zhaires
                                DAmpy = Ampy_RM-Ampy_zhaires
                                DAmpz = Ampz_RM-Ampz_zhaires
                                DAmpx_rel = (Ampx_RM-Ampx_zhaires)/Ampx_zhaires
                                DAmpy_rel = (Ampy_RM-Ampy_zhaires)/Ampy_zhaires
                                DAmpz_rel = (Ampz_RM-Ampz_zhaires)/Ampz_zhaires

                                DAmpx_rel = np.reshape(DAmpx_rel,np.size(DAmpx_rel))
                                DAmpy_rel = np.reshape(DAmpy_rel,np.size(DAmpy_rel))
                                DAmpz_rel = np.reshape(DAmpz_rel,np.size(DAmpz_rel))
                                DAmpx = np.reshape(DAmpx,np.size(DAmpx))
                                DAmpy = np.reshape(DAmpy,np.size(DAmpy))
                                DAmpz = np.reshape(DAmpz,np.size(DAmpz))

                                ind = zen_ind[np.where(showerID_sh==showerID)][0]
                                DAmpx_rel_all[ind] = np.concatenate((DAmpx_rel_all[ind],DAmpx_rel))
                                DAmpy_rel_all[ind] = np.concatenate((DAmpy_rel_all[ind],DAmpy_rel))
                                DAmpz_rel_all[ind] = np.concatenate((DAmpz_rel_all[ind],DAmpz_rel))

                                DAmpx_all[ind] = np.concatenate((DAmpx_all[ind],DAmpx))
                                DAmpy_all[ind] = np.concatenate((DAmpy_all[ind],DAmpy))
                                DAmpz_all[ind] = np.concatenate((DAmpz_all[ind],DAmpz))
                            else:
                                print('no antenna')
                                pass

#    DAmpx_rel_all = [x for x in DAmpx_rel_all if not np.isnan(x)]
#    DAmpy_rel_all = [x for x in DAmpy_rel_all if not np.isnan(x)]
#    DAmpz_rel_all = [x for x in DAmpz_rel_all if not np.isnan(x)]
#    DAmpx_all = [x for x in DAmpx_all if not np.isnan(x)]
#    DAmpy_all = [x for x in DAmpy_all if not np.isnan(x)]
#    DAmpz_all = [x for x in DAmpz_all if not np.isnan(x)]


    DAmpx_rel_min = [np.amin(DAmpx_rel_all[ix]) if len(DAmpx_rel_all[ix])!=0 else 0. for ix in range(0,len(DAmpx_rel_all))]
    DAmpx_rel_max = [np.amax(DAmpx_rel_all[ix]) if len(DAmpx_rel_all[ix])!=0 else 1. for ix in range(0,len(DAmpx_rel_all))]
    DAmpy_rel_min = [np.amin(DAmpy_rel_all[ix]) if len(DAmpy_rel_all[ix])!=0 else 0. for ix in range(0,len(DAmpy_rel_all))]
    DAmpy_rel_max = [np.amax(DAmpy_rel_all[ix]) if len(DAmpy_rel_all[ix])!=0 else 1. for ix in range(0,len(DAmpy_rel_all))]
    DAmpz_rel_min = [np.amin(DAmpz_rel_all[ix]) if len(DAmpz_rel_all[ix])!=0 else 0. for ix in range(0,len(DAmpz_rel_all))]
    DAmpz_rel_max = [np.amax(DAmpz_rel_all[ix]) if len(DAmpz_rel_all[ix])!=0 else 1. for ix in range(0,len(DAmpz_rel_all))]
    DAmpx_min = [np.amin(DAmpx_all[ix]) if len(DAmpx_all[ix])!=0 else -1. for ix in range(0,len(DAmpx_all))]
    DAmpx_max = [np.amax(DAmpx_all[ix]) if len(DAmpx_all[ix])!=0 else 1. for ix in range(0,len(DAmpx_all))]
    DAmpy_min = [np.amin(DAmpy_all[ix]) if len(DAmpy_all[ix])!=0 else -1. for ix in range(0,len(DAmpy_all))]
    DAmpy_max = [np.amax(DAmpy_all[ix]) if len(DAmpy_all[ix])!=0 else 1. for ix in range(0,len(DAmpy_all))]
    DAmpz_min = [np.amin(DAmpz_all[ix]) if len(DAmpz_all[ix])!=0 else -1. for ix in range(0,len(DAmpz_all))]
    DAmpz_max = [np.amax(DAmpz_all[ix]) if len(DAmpz_all[ix])!=0 else 1. for ix in range(0,len(DAmpz_all))]

#    DAmpx_rel_bin = np.arange(DAmpx_rel_min,DAmpx_rel_max,(DAmpx_rel_max-DAmpx_rel_min)/1000.)
#    DAmpy_rel_bin = np.arange(DAmpy_rel_min,DAmpy_rel_max,(DAmpy_rel_max-DAmpy_rel_min)/1000.)
#    DAmpz_rel_bin = np.arange(DAmpz_rel_min,DAmpz_rel_max,(DAmpz_rel_max-DAmpz_rel_min)/1000.)
#    DAmpx_bin = np.arange(DAmpx_min,DAmpx_max,(DAmpx_max-DAmpx_min)/1000.)
#    DAmpy_bin = np.arange(DAmpy_min,DAmpy_max,(DAmpy_max-DAmpy_min)/1000.)
#    DAmpz_bin = np.arange(DAmpz_min,DAmpz_max,(DAmpz_max-DAmpz_min)/1000.)

    nbin=20 #'auto' #30
    fig1 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    fracwidth = np.floor((np.amax(DAmpx_rel_max)-np.amin(DAmpx_rel_min))/((nbin+1)*(len(DAmpx_rel_all)+1))*1e3)/1e3
    hist_DAmpx_rel_all,bin_DAmpx_rel_all,tmp = pl.hist(DAmpx_rel_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpx_rel_bin,hist_fit_Dampx_rel,'--r',label=label_rel_x,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{x}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    fracwidth = np.floor((np.amax(DAmpy_rel_max)-np.amin(DAmpy_rel_min))/((nbin+1)*(len(DAmpy_rel_all)+1))*1e3)/1e3
    hist_DAmpy_rel_all,bin_DAmpy_rel_all,tmp = pl.hist(DAmpy_rel_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpy_rel_bin,hist_fit_Dampy_rel,'--r',label=label_rel_y,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{y}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    fracwidth = np.floor((np.amax(DAmpz_rel_max)-np.amin(DAmpz_rel_min))/((nbin+1)*(len(DAmpz_rel_all)+1))*1e3)/1e3
    hist_DAmpz_rel_all,bin_DAmpz_rel_all,tmp = pl.hist(DAmpz_rel_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpz_rel_bin,hist_fit_Dampz_rel,'--r',label=label_rel_z,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Relative difference on E$_{z}$')
    pl.ylabel('Rate')
    #pl.subplot(224)
    #pl.hist([DAmpx_rel_all,DAmpy_rel_all,DAmpz_rel_all],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    #pl.xlabel('Relative difference on E')
    #pl.ylabel('Rate')
    #pl.legend()
    figname = figdir+'/RelDeltaE_high.png'
    pl.savefig(figname,dpi=350)

    fig2 = pl.figure(figsize=(21.0,9.7)) #(5*3.13,3.1*3.13)) #figsize=(5*3.13,3.9*3.13))
    pl.subplot(221)
    fracwidth = np.floor((np.amax(DAmpx_max)-np.amin(DAmpx_min))/((nbin+1)*(len(DAmpx_all)+1))*1e3)/1e3
    hist_DAmpx_all,bin_DAmpx_all,tmp = pl.hist(DAmpx_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpx_bin,hist_fit_Dampx,'--r',label=labelx,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{x}$')
    pl.ylabel('Rate')
    pl.subplot(222)
    fracwidth = np.floor((np.amax(DAmpy_max)-np.amin(DAmpy_min))/((nbin+1)*(len(DAmpy_all)+1))*1e3)/1e3
    hist_DAmpy_all,bin_DAmpy_all,tmp = pl.hist(DAmpy_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpy_bin,hist_fit_Dampy,'--r',label=labely,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{y}$')
    pl.ylabel('Rate')
    pl.subplot(223)
    fracwidth = np.floor((np.amax(DAmpz_max)-np.amin(DAmpz_min))/((nbin+1)*(len(DAmpz_all)+1))*1e3)/1e3
    hist_DAmpz_all,bin_DAmpz_all,tmp = pl.hist(DAmpz_all,bins=nbin,normed=1,width=fracwidth)
    try:
        pl.plot(DAmpz_bin,hist_fit_Dampz,'--r',label=labelz,linewidth=2)
        pl.legend()
    except:
        pass
    pl.xlabel('Absolute difference on E$_{z}$')
    pl.ylabel('Rate')
    #pl.subplot(224)
    #pl.hist([DAmpx_all,DAmpy_all,DAmpz_all],bins=nbin,normed=1,label=[r'$\rm \Delta E_{x}$',r'$\rm \Delta E_{y}$',r'$\rm \Delta E_{z}$'])
    #pl.xlabel('Absolute difference on E')
    #pl.ylabel('Rate')
    #pl.legend()
    #figname = figdir+'/AbsDeltaE.png'
    #pl.savefig(figname,dpi=350)

    '''
    print('\n mu Rel DEx',np.round(np.mean(DAmpx_rel_all),2),'mu Rel DEy',np.round(np.mean(DAmpy_rel_all),2),'mu Rel DEz',np.round(np.mean(DAmpz_rel_all),2))
    print(' std Rel DEx',np.round(np.std(DAmpx_rel_all),2),'std Rel DEy',np.round(np.std(DAmpy_rel_all),2),'std Rel DEz',np.round(np.std(DAmpz_rel_all),2))
    print(' med Rel DEx',np.round(np.median(DAmpx_rel_all),2),'med Rel DEy',np.round(np.median(DAmpy_rel_all),2),'med Rel DEz',np.round(np.median(DAmpz_rel_all),2))
    print(' pos max Rel DEx',np.round(bin_DAmpx_rel_all[np.argmax(hist_DAmpx_rel_all)],2),'pos max Rel DEy',np.round(bin_DAmpy_rel_all[np.argmax(hist_DAmpy_rel_all)],2),'pos max Rel DEz',np.round(bin_DAmpz_rel_all[np.argmax(hist_DAmpz_rel_all)],2))
    print('\n mu Abso DEx',np.round(np.mean(DAmpx_all),2),'mu Abso DEy',np.round(np.mean(DAmpy_all),2),'mu Abso DEz',np.round(np.mean(DAmpz_all),2))
    print(' std Abso DEx',np.round(np.std(DAmpx_all),2),'std Abso DEy',np.round(np.std(DAmpy_all),2),'std Abso DEz',np.round(np.std(DAmpz_all),2))
    print(' med Abso DEx',np.round(np.median(DAmpx_all),2),'med Abso DEy',np.round(np.median(DAmpy_all),2),'med Abso DEz',np.round(np.median(DAmpz_all),2))
    print(' pos max Abso DEx',np.round(bin_DAmpx_all[np.argmax(hist_DAmpx_all)],2),'pos max Abso DEy',np.round(bin_DAmpy_all[np.argmax(hist_DAmpy_all)],2),'pos max Abso DEz',np.round(bin_DAmpz_all[np.argmax(hist_DAmpz_all)],2))
    print(' ')
    '''
    if DISPLAY:
        pl.show()
    pl.close(fig1)
    pl.close(fig2)

    return

############################################################################################################
############################################################################################################
############################################################################################################
if __name__ == '__main__':

    ROOT_RM = "/Users/nrenault/Desktop/GRAND/RadioMorphing/simus_RM/"
    ROOT_zhaires = "/Users/nrenault/Desktop/GRAND/RadioMorphing/simus_zhaires/"

    Nserie=[1] #[1,2]
    EE='1E+10' #Neutrino energy e in GeV
    hz=3000 #Mountain Height
    sep=[500] #[500,500,500,250,250,250,250] #separation between antennas

    #Dd=[20000,30000,40000,60000,80000,100000] #distance to decay point in m
    #alfa=[10,15,20,45,90] #slope angle in degrees
    #Etarget = np.array([5E+17, 1E+18, 1E+19],dtype=str) #eV #1e19 #1e18 # 5e17
    Dd=[40000] #distance to decay point in m
    alfa=[10] #slope angle in degrees
    Etarget = np.array([1E+19],dtype=str) #eV #1e19 #1e18 # 5e17

    fct_sel = int(sys.argv[1])
    if fct_sel==1:
        histo_DEpp(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires)
    elif fct_sel==2:
        histo_DEpp_threshold(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires)
    elif fct_sel==3:
        histo_Ntrig(Dd,alfa,hz,Etarget,sep,ROOT_RM,ROOT_zhaires)
    elif fct_sel==4:
        histo_DEpp_EnAzZen(Dd,alfa,hz,EE,sep,Etarget,Nserie,ROOT_RM,ROOT_zhaires)
    else:
        print("Set the function number you want to run : 1,2 or 3")
