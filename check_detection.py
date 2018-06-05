#import matplotlib
#matplotlib.use('Agg')

import os, glob
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys
import numpy as np
import linecache
import StringIO

###########################
###    Test arguments   ###
###########################
if (len(sys.argv)<1 or len(sys.argv)>4):
    print """\
        This script will read the voltage for all antennas, compute the p2p voltage, and test how many are above thresholdself.
        A DISPLAY variable can be set to one to produce the voltage and trigger maps.
        The script writes the results in a txt file.

        Usage:  python check_detection.py [toymodel, RM, CR or GP300] [opt: low freq] [opt: high freq]
        """
    sys.exit(1)

threshold_cons=150. #noise=15muV => case with 10*noise
threshold_aggr=50. #noise=15muV => case with 10/3*noise

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

####################################################################################
def read_voltage(trace,suffix):
    output_folder=os.path.dirname(trace)
    tmp = os.path.basename(trace)
    tmp1 = tmp.replace('/','.').split('.')[-2][1:]

    volt_file=output_folder+'/out_'+tmp1+suffix+'.txt'
    if os.path.exists(volt_file):
        s = open(volt_file).read().replace('  ',' ').replace('\t',' ')
        time,Vns,Vew,Vup = np.loadtxt(StringIO.StringIO(s),delimiter=' ',usecols=[0,1,2,3],unpack=True)
        time = np.array(time,dtype=float)
        Vew = np.array(Vew,dtype=float)
        Vns = np.array(Vns,dtype=float)
        Vup = np.array(Vup,dtype=float)
        #Vp2p = [np.max(Vew)+np.abs(np.min(Vew)),np.max(Vns)+np.abs(np.min(Vns)),np.max(Vup)+np.abs(np.min(Vup)),np.max(Vew+Vns+Vup)+np.abs(np.min(Vew+Vns+Vup))]
        Vp2p = [np.max(Vew)+np.abs(np.min(Vew)),np.max(Vns)+np.abs(np.min(Vns)),np.max(Vup)+np.abs(np.min(Vup)),np.max(Vew+Vns)+np.abs(np.min(Vew+Vns))]
    else:
        Vp2p = [0.,0.,0.,0.]

    tmp1 = tmp.replace('/','.').split('.')[0]
    numberline = int(tmp1.split('a')[1]) + 1
    line = linecache.getline(output_folder+'/antpos.dat', numberline)
    antpos = map(float, line.split())
    return np.concatenate((Vp2p,antpos))

####################################################################################
def compute_toymodel(Dd,alfa,hz,Nserie,EE,sep,Etarget,ROOT,suffix):
    for iserie in Nserie:
        for ie in range(len(Etarget)):
            ratioE = float(Etarget[ie])/1E+19

            #trig_file=ROOT+'detection/detection_count.txt'
            #shower_file_tot=ROOT+'detection/detection_count_tot_EE'+Etarget[ie]+suffix+'.txt'
            data_trig = []
            #data_config_tot=[]
            for iDd in Dd:
                for ialfa in range(np.size(alfa)):
                    if alfa[ialfa]>0 or (alfa[ialfa]==0 and iDd==20000):
                        data_config=[]
                        test_folder=ROOT+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+'/'
                        config_folder=ROOT+'detection/'+'EE'+Etarget[ie]+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+'/'
                        shower_file=ROOT+'detection/detection_count_EE'+Etarget[ie]+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'_serie'+str(iserie)+suffix+'.txt'
                        showers = glob.glob(test_folder+'[0-9][0-9]*[0-9]/')
                        #showers=showers[103:104]

                        #Start loop on shower
                        for shower in showers:
                            showerID = shower.replace('/','.').split('.')[-2]
                            trace_files=glob.glob(shower+'split/*.trace')
                            Nant=len(trace_files)

                            #Start loop on antennas
                            if Nant>0:
                                Vp2p_antpos = np.array(map(lambda x:read_voltage(x,suffix),trace_files))
                                Vp2p_ew = Vp2p_antpos[:,0]*ratioE
                                Vp2p_ns = Vp2p_antpos[:,1]*ratioE
                                Vp2p_up = Vp2p_antpos[:,2]*ratioE
                                Vp2p_tot = Vp2p_antpos[:,3]*ratioE
                                antpos = Vp2p_antpos[:,4:9]

                                ind_ew_cons = np.where(Vp2p_ew>=threshold_cons)
                                ind_ns_cons = np.where(Vp2p_ns>=threshold_cons)
                                ind_up_cons = np.where(Vp2p_up>=threshold_cons)
                                ind_tot_cons = np.where(Vp2p_tot>=threshold_cons)
                                ind_ew_aggr = np.where(Vp2p_ew>=threshold_aggr)
                                ind_ns_aggr = np.where(Vp2p_ns>=threshold_aggr)
                                ind_up_aggr = np.where(Vp2p_up>=threshold_aggr)
                                ind_tot_aggr = np.where(Vp2p_tot>=threshold_aggr)
                                Nant_ew_cons = np.size(ind_ew_cons)
                                Nant_ns_cons = np.size(ind_ns_cons)
                                Nant_up_cons = np.size(ind_up_cons)
                                Nant_tot_cons = np.size(ind_tot_cons)
                                Nant_ew_aggr = np.size(ind_ew_aggr)
                                Nant_ns_aggr = np.size(ind_ns_aggr)
                                Nant_up_aggr = np.size(ind_up_aggr)
                                Nant_tot_aggr = np.size(ind_tot_aggr)

                                data_config.append((iDd/1e3,alfa[ialfa],showerID,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))
                                data_config_tot.append((iDd/1e3,alfa[ialfa],showerID,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))

                                #Plot
                                DISPLAY=0
                                if DISPLAY:
                                    if not os.path.exists(config_folder):
                                        os.makedirs(config_folder)
                                    plot_maps(Vp2p_ew,Vp2p_ns,Vp2p_up,antpos,config_folder,showerID,suffix)
                            else:
                                data_config.append((iDd/1e3,alfa[ialfa],showerID,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                                #data_config_tot.append((iDd/1e3,alfa[ialfa],showerID,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                            #End loop on antennas
                        #End loop on showers

                        data_config = np.array(data_config,dtype=str)
                        np.savetxt(shower_file,data_config,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s   %s   %s ', header="Distance [km]   Slope [deg]   ShowerId   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
            #np.savetxt(shower_file_tot,data_config_tot,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s   %s   %s     %s ', header="Distance [km]   Slope [deg]   ShowerId   New_conservative   Nns_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_tot) ")
    return

####################################################################################
def compute_RM(Dd,alfa,hz,EE,sep,Etarget,ROOT):
    for iserie in Nserie:
        for ie in range(len(Etarget)):
            ratioE = float(Etarget[ie])/1E+19

            data_trig = []
            data_config_tot=[]
            for iDd in Dd:
                for ialfa in range(np.size(alfa)):
                    if alfa[ialfa]>0 or (alfa[ialfa]==0 and iDd==20000):
                        data_config=[]
                        test_folder=ROOT+'EE'+EE+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'/'
                        config_folder=ROOT+'detection/'+'EE'+Etarget[ie]+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'/'
                        shower_file=ROOT+'detection/detection_count_EE'+Etarget[ie]+'_D'+str(iDd)+'_alpha'+str(alfa[ialfa])+'_height'+str(hz)+'_sep'+str(sep[ialfa])+'.txt'
                        showers = glob.glob(test_folder+'[0-9][0-9]*[0-9]/')

                        #Start loop on shower
                        for shower in showers:
                            showerID = shower.replace('/','.').split('.')[-2]
                            trace_files=glob.glob(shower+'*.trace')
                            Nant=len(trace_files)

                            #Start loop on antennas
                            if Nant>0:
                                Vp2p_antpos = np.array(map(lambda x:read_voltage(x),trace_files))
                                Vp2p_ew = Vp2p_antpos[:,0]*ratioE
                                Vp2p_ns = Vp2p_antpos[:,1]*ratioE
                                Vp2p_up = Vp2p_antpos[:,2]*ratioE
                                Vp2p_tot = Vp2p_antpos[:,3]*ratioE
                                antpos = Vp2p_antpos[:,4:9]

                                ind_ew_cons = np.where(Vp2p_ew>=threshold_cons)
                                ind_ns_cons = np.where(Vp2p_ns>=threshold_cons)
                                ind_up_cons = np.where(Vp2p_up>=threshold_cons)
                                ind_tot_cons = np.where(Vp2p_tot>=threshold_cons)
                                ind_ew_aggr = np.where(Vp2p_ew>=threshold_aggr)
                                ind_ns_aggr = np.where(Vp2p_ns>=threshold_aggr)
                                ind_up_aggr = np.where(Vp2p_up>=threshold_aggr)
                                ind_tot_aggr = np.where(Vp2p_tot>=threshold_aggr)
                                Nant_ew_cons = np.size(ind_ew_cons)
                                Nant_ns_cons = np.size(ind_ns_cons)
                                Nant_up_cons = np.size(ind_up_cons)
                                Nant_tot_cons = np.size(ind_tot_cons)
                                Nant_ew_aggr = np.size(ind_ew_aggr)
                                Nant_ns_aggr = np.size(ind_ns_aggr)
                                Nant_up_aggr = np.size(ind_up_aggr)
                                Nant_tot_aggr = np.size(ind_tot_aggr)

                                data_config.append((iDd/1e3,ialfa,showerID,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))
                                data_config_tot.append((iDd/1e3,ialfa,showerID,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))

                                #Plot
                                DISPLAY=0
                                if DISPLAY:
                                    if not os.path.exists(config_folder):
                                        os.makedirs(config_folder)
                                    plot_maps(Vp2p_ew,Vp2p_ns,Vp2p_up,antpos,config_folder,showerID)

                            else:
                                data_config.append((iDd/1e3,ialfa,showerID,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                                data_config_tot.append((iDd/1e3,ialfa,showerID,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                            #End loop on antennas
                        #End loop on showers

                        data_config = np.array(data_config,dtype=str)
                        np.savetxt(shower_file,data_config,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s   %s   %s ', header="Distance [km]   Slope [deg]   ShowerId   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
    return

##########################################################################################################
def compute_CR(az,zen,Etarget,Nevt,ROOT,suffix):
    shower_file_tot=ROOT+'detection/detection_count'+suffix+'.txt'
    data_config_tot=[]
    for ievt in Nevt:
        for ie in range(len(Etarget)):
            ratioE = 1. #float(Etarget[ie])/1E+19
            for iaz in az:
                for izen in zen:
                    data_config=[]
                    task = 'EE'+Etarget[ie]+'-az'+str(iaz)+'-zen'+str(izen)+'_evt'+str(ievt)
                    test_folder=ROOT+'/evt'+str(ievt)+'/'+task+'/'
                    config_folder=ROOT+'detection/evt'+str(ievt)+'/'+task+'/'
                    shower_file=ROOT+'detection/detection_count_'+task+suffix+'.txt'
                    trace_files=glob.glob(test_folder+'split/*.trace')
                    Nant=len(trace_files)

                    #Start loop on antennas
                    if Nant>0:
                        Vp2p_antpos = np.array(map(lambda x:read_voltage(x,suffix),trace_files))
                        Vp2p_ew = Vp2p_antpos[:,0]*ratioE
                        Vp2p_ns = Vp2p_antpos[:,1]*ratioE
                        Vp2p_up = Vp2p_antpos[:,2]*ratioE
                        Vp2p_tot = Vp2p_antpos[:,3]*ratioE
                        antpos = Vp2p_antpos[:,4:9]

                        ind_ew_cons = np.where(Vp2p_ew>=threshold_cons)
                        ind_ns_cons = np.where(Vp2p_ns>=threshold_cons)
                        ind_up_cons = np.where(Vp2p_up>=threshold_cons)
                        ind_tot_cons = np.where(Vp2p_tot>=threshold_cons)
                        ind_ew_aggr = np.where(Vp2p_ew>=threshold_aggr)
                        ind_ns_aggr = np.where(Vp2p_ns>=threshold_aggr)
                        ind_up_aggr = np.where(Vp2p_up>=threshold_aggr)
                        ind_tot_aggr = np.where(Vp2p_tot>=threshold_aggr)
                        Nant_ew_cons = np.size(ind_ew_cons)
                        Nant_ns_cons = np.size(ind_ns_cons)
                        Nant_up_cons = np.size(ind_up_cons)
                        Nant_tot_cons = np.size(ind_tot_cons)
                        Nant_ew_aggr = np.size(ind_ew_aggr)
                        Nant_ns_aggr = np.size(ind_ns_aggr)
                        Nant_up_aggr = np.size(ind_up_aggr)
                        Nant_tot_aggr = np.size(ind_tot_aggr)

                        data_config.append((Etarget[ie],iaz,izen,ievt,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))
                        data_config_tot.append((Etarget[ie],iaz,izen,ievt,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))

                        #Plot
                        DISPLAY=0
                        if DISPLAY:
                            if not os.path.exists(config_folder):
                                os.makedirs(config_folder)
                            plot_maps(Vp2p_ew,Vp2p_ns,Vp2p_up,antpos,config_folder,task,suffix)

                    else:
                        data_config.append((Etarget[ie],iaz,izen,ievt,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                        data_config_tot.append((Etarget[ie],iaz,izen,ievt,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                    #End loop on antennas

                    data_config = np.array(data_config,dtype=str)
                    #np.savetxt(shower_file,data_config,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s   %s   %s ', header="Energy [eV]    Azimuth_GRAND [deg]   Zenith_GRAND [deg]   EventNumber   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
    data_config_tot = np.array(data_config_tot,dtype=str)
    np.savetxt(shower_file_tot,data_config_tot,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s   %s   %s ', header="Energy [eV]    Azimuth_GRAND [deg]   Zenith_GRAND [deg]   EventNumber   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
    return

##########################################################################################################
def compute_GP300(Etarget,Nevt,ROOT,suffix,DISPLAY):
    det_dir = ROOT+'detection/'
    if not os.path.exists(det_dir):
        os.makedirs(det_dir)
    shower_file_tot=ROOT+'detection/detection_count'+suffix+'.txt'
    data_config_tot=[]
    for ie in range(len(Etarget)):
        for ievt in Nevt:
            ratioE = 1. #10**float(Etarget[ie])/1E+19
            data_config=[]
            task = 'EE1E'+Etarget[ie]+'_evt'+str(ievt)
            test_folder=ROOT+'/'+task+'/'
            config_folder=ROOT+'detection/'+task+'/'
            shower_file=ROOT+'detection/detection_count_'+task+suffix+'.txt'
            trace_files=glob.glob(test_folder+'split/*.trace')
            Nant=len(trace_files)

            #Start loop on antennas
            if Nant>0:
                Vp2p_antpos = np.array(map(lambda x:read_voltage(x,suffix),trace_files))
                Vp2p_ew = Vp2p_antpos[:,0]*ratioE
                Vp2p_ns = Vp2p_antpos[:,1]*ratioE
                Vp2p_up = Vp2p_antpos[:,2]*ratioE
                Vp2p_tot = Vp2p_antpos[:,3]*ratioE
                antpos = Vp2p_antpos[:,4:9]

                ind_ew_cons = np.where(Vp2p_ew>=threshold_cons)
                ind_ns_cons = np.where(Vp2p_ns>=threshold_cons)
                ind_up_cons = np.where(Vp2p_up>=threshold_cons)
                ind_tot_cons = np.where(Vp2p_tot>=threshold_cons)
                ind_ew_aggr = np.where(Vp2p_ew>=threshold_aggr)
                ind_ns_aggr = np.where(Vp2p_ns>=threshold_aggr)
                ind_up_aggr = np.where(Vp2p_up>=threshold_aggr)
                ind_tot_aggr = np.where(Vp2p_tot>=threshold_aggr)
                Nant_ew_cons = np.size(ind_ew_cons)
                Nant_ns_cons = np.size(ind_ns_cons)
                Nant_up_cons = np.size(ind_up_cons)
                Nant_tot_cons = np.size(ind_tot_cons)
                Nant_ew_aggr = np.size(ind_ew_aggr)
                Nant_ns_aggr = np.size(ind_ns_aggr)
                Nant_up_aggr = np.size(ind_up_aggr)
                Nant_tot_aggr = np.size(ind_tot_aggr)

                data_config.append((Etarget[ie],ievt,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))
                data_config_tot.append((Etarget[ie],ievt,Nant_ew_cons,Nant_ns_cons,Nant_up_cons,Nant_tot_cons,Nant_ew_aggr,Nant_ns_aggr,Nant_up_aggr,Nant_tot_aggr,np.max(Vp2p_ew),np.max(Vp2p_ns),np.max(Vp2p_up),np.max(Vp2p_tot)))

                #Plot
                if DISPLAY:
                    if not os.path.exists(config_folder):
                        os.makedirs(config_folder)
                    plot_maps(Vp2p_ew,Vp2p_ns,Vp2p_up,Vp2p_tot,antpos,config_folder,task,suffix)

            else:
                data_config.append((Etarget[ie],ievt,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
                data_config_tot.append((Etarget[ie],ievt,0,0,0,0,0,0,0,0,0.,0.,0.,0.))
            #End loop on antennas

            data_config = np.array(data_config,dtype=str)
            #np.savetxt(shower_file,data_config,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s   %s   %s ', header="Energy [eV]    Azimuth_GRAND [deg]   Zenith_GRAND [deg]   EventNumber   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
    data_config_tot = np.array(data_config_tot,dtype=str)
    np.savetxt(shower_file_tot,data_config_tot,fmt='%s     %s     %s     %s     %s     %s     %s     %s     %s     %s     %s    %s   %s   %s ', header="Energy [eV]   EventNumber   New_conservative   Nns_conservative   Nup_conservative   Ntot_conservative   New_aggressive   Nns_aggressive   Nup_aggressive   Ntot_aggressive   max(Vp2p_ew)   max(Vp2p_ns)   max(Vp2p_up)   max(Vp2p_tot) ")
    return

##########################################################################################################
def plot_maps(Vp2p_ew,Vp2p_ns,Vp2p_up,Vp2p_tot,antpos,config_folder,showerID,suffix):
    binmap_ew_cons = np.zeros(np.shape(Vp2p_ew))
    binmap_ew_cons[Vp2p_ew>=threshold_cons]=1
    binmap_ew_aggr = np.zeros(np.shape(Vp2p_ew))
    binmap_ew_aggr[Vp2p_ew>=threshold_aggr]=1

    binmap_ns_cons = np.zeros(np.shape(Vp2p_ns))
    binmap_ns_cons[Vp2p_ns>=threshold_cons]=1
    binmap_ns_aggr = np.zeros(np.shape(Vp2p_ns))
    binmap_ns_aggr[Vp2p_ns>=threshold_aggr]=1

    binmap_up_cons = np.zeros(np.shape(Vp2p_up))
    binmap_up_cons[Vp2p_up>=threshold_cons]=1
    binmap_up_aggr = np.zeros(np.shape(Vp2p_up))
    binmap_up_aggr[Vp2p_up>=threshold_aggr]=1

    binmap_tot_cons = np.zeros(np.shape(Vp2p_tot))
    binmap_tot_cons[Vp2p_tot>=threshold_cons]=1
    binmap_tot_aggr = np.zeros(np.shape(Vp2p_tot))
    binmap_tot_aggr[Vp2p_tot>=threshold_aggr]=1

    binmap = ListedColormap(['white', 'black'], 'indexed')
    dar=(np.max(antpos[:,0])-np.min(antpos[:,0]))/(np.max(antpos[:,1])-np.min(antpos[:,1]))/7
    if dar==0:
        dar=1

    #EW
    fig1 = pl.figure(1)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=Vp2p_ew,cmap='jet',label='Vew',edgecolors='none')
    pl.title('Voltage EW [muV]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    pl.colorbar()
    #adjustFigAspect(fig1,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_voltmap_ew'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig1)

    fig2 = pl.figure(2)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_ew_cons,cmap=binmap,label='Vew',vmin=0.,vmax=1.)
    pl.title('Triggered EW (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig2,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_cons_ew'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig2)

    fig3 = pl.figure(3)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_ew_aggr,cmap=binmap,label='Vew',vmin=0.,vmax=1.)
    pl.title('Triggered EW (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig3,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_aggr_ew'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig3)

    #NS
    fig4 = pl.figure(4)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=Vp2p_ns,cmap='jet',label='Vns',edgecolors='none')
    pl.title('Voltage NS [muV]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    pl.colorbar()
    #adjustFigAspect(fig4,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_voltmap_ns'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig4)

    fig5 = pl.figure(5)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_ns_cons,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered NS (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig5,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_cons_ns'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig5)

    fig6 = pl.figure(6)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_ns_aggr,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered NS (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig6,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_aggr_ns'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig6)

    #UP
    fig7 = pl.figure(7)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=Vp2p_up,cmap='jet',label='Vns',edgecolors='none')
    pl.title('Voltage Vertical Arm [muV]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    pl.colorbar()
    #adjustFigAspect(fig4,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_voltmap_up'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig7)

    fig8 = pl.figure(8)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_up_cons,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered Vertical Arm (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig5,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_cons_up'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig8)

    fig9 = pl.figure(9)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_up_aggr,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered Vertical Arm (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig6,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_aggr_up'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig9)

    #UP
    fig10 = pl.figure(10)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=Vp2p_tot,cmap='jet',label='Vns',edgecolors='none')
    pl.title('Total Voltage [muV]')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    pl.colorbar()
    #adjustFigAspect(fig4,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_voltmap_tot'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig10)

    fig11 = pl.figure(11)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_tot_cons,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered Total Voltage (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig5,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_cons_tot'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig11)

    fig12 = pl.figure(12)
    xlbl='X [m]'
    ylbl='Y [m]'
    pl.scatter(antpos[:,0],antpos[:,1],c=binmap_tot_aggr,cmap=binmap,label='Vns',vmin=0.,vmax=1.)
    pl.title('Triggered Total Voltage (True/False)')
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    #adjustFigAspect(fig6,aspect=dar)
    pl.axis('equal')
    figname = config_folder+'/'+showerID+'_trigmap_aggr_tot'+suffix+'.png' #_50-200MHz
    pl.savefig(figname,dpi=500)
    pl.close(fig12)
    return

##########################################################################################################
if __name__ == '__main__':
    analysis = sys.argv[1]

    if analysis=='toymodel':
        ROOT = "/Users/nrenault/Desktop/GRAND/RadioMorphing/simus_zhaires/"
        #ROOT = "/data75/renault/ToyModel/"
        #ROOT = "/Users/nrenault/Desktop/GRAND/ToyModel/"

        Nserie=[1] #[1,2]
        EE='1E+10' #Neutrino energy e in GeV
        hz=3000 #Mountain Height
        sep=[500,500,500,250,250,250,250] #separation between antennas

        Dd=[20000,30000,40000,60000,80000,100000] #distance to decay point in m
        alfa=[10,15,20,45,90] #slope angle in degrees
        Etarget = np.array([5E+17, 1E+18, 1E+19],dtype=str) #eV #1e19 #1e18 # 5e17

        '''
        try:
            lowcut = int(sys.argv[1])
            highcut = int(sys.argv[2])
            suffix = '_'+str(lowcut)+'-'+str(highcut)+'MHz'
        except:
          lowcut=0
          highcut=0
          suffix=''
        '''

        suffix = '_50-200MHz'
        compute_toymodel(Dd,alfa,hz,Nserie,EE,sep,Etarget,ROOT,suffix)

    if analysis=='RM':
        ROOT = "/Users/nrenault/Desktop/GRAND/RadioMorphing/simus_RM/"
        #ROOT = "/data75/renault/ToyModel/"

        Nserie=[1] #[1,2]
        EE='1E+10' #Neutrino energy e in GeV
        hz=3000 #Mountain Height
        sep=[500] #[500,500,500,250,250,250,250] #separation between antennas

        Dd=[20000,30000,40000,60000,80000,100000] #distance to decay point in m
        alfa=[10,15,20,45,90] #slope angle in degrees
        Etarget = np.array([5E+17, 1E+18, 1E+19],dtype=str) #eV #1e19 #1e18 # 5e17

        compute_RM(Dd,alfa,hz,EE,sep,Etarget,ROOT)

    elif analysis=='CR':
        ROOT = "/Users/nrenault/Desktop/GRAND/CR_v1/simus_test/"
        #ROOT = "/data75/renault/CR/simus/"
        az=[0,45,90,135,180] #distance to decay point in m
        zen=[95,100,105,110,115,120] #slope angle in degrees
        Etarget = np.array(['1E17.5', '1E18.0', '1E18.5', '1E19.0', '1E19.5'],dtype=str) #eV #1e19 #1e18 # 5e17
        Nevt = range(0,10)

        try:
          lowcut = int(sys.argv[1])
          highcut = int(sys.argv[2])
          suffix = '_'+str(lowcut)+'-'+str(highcut)+'MHz'
        except:
          lowcut=0
          highcut=0
          suffix=''

        compute_CR(az,zen,Etarget,Nevt,ROOT,suffix)

    elif analysis=='GP300':
        DISPLAY=0
        ROOT = "/Users/nrenault/Desktop/GRAND/CRs_GP300/simus_set3/"
        #ROOT = "/data75/renault/CRs_GP300/simus/voltage_computed/"

        Etarget = np.array(['17.0','17.5','18.0', '18.5','19.0', '19.5'],dtype=str)# , '19.0', '19.5'],dtype=str) #eV #1e19 #1e18 # 5e17
        #Etarget = np.array(['16.5'],dtype=str)# , '19.0', '19.5'],dtype=str) #eV #1e19 #1e18 # 5e17
        Nevt = range(1,301)

        '''
        try:
          lowcut = int(sys.argv[1])
          highcut = int(sys.argv[2])
          suffix = '_'+str(lowcut)+'-'+str(highcut)+'MHz'
        except:
          lowcut=0
          highcut=0
          suffix=''
        '''

        suffix = '_50-200MHz'
        compute_GP300(Etarget,Nevt,ROOT,suffix,DISPLAY)

    print '\n'
    print 'Finished'
