''' this script bases on a script provided by S.LeCoz to filter, mimic digitization and add noise after applying an antenna response
    It can be only applied for complete array or single antenna positions having array.dat and antpos.dat as file
    
    hand over path to folder containing somulations and antenna files as argument
    -> out_*.txt as output from computevoltage.py
    -> antpos.dat containing all antenna position with ID of all antennas which where included in the simulations
    -> array.dat containing the antenna position of the whole field
    
    => output: fake_*.dat with the antenna ID following the lines in array.dat : time in s, EW/nuV, NS/muV (to be checked)
    if antenna position in array.dat was not icnluded in the simulations, this script creates an empty trace and does the filtering, Digitization and the noise adding with it
    
    => prepare input for NEURAL NETWORK
''' 

import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from random import *
import matplotlib.pyplot as plt
import sys


def Filtering(v,tstep,FREQMIN,FREQMAX):
    F=rfftfreq(len(v))/tstep #len(v) points between 0 and 1/2tstep
    V=rfft(v)
    V[F<FREQMIN]=0
    V[F>FREQMAX]=0
    return irfft(V)

def Digitization(v,t,tstep,TSAMPLING,SAMPLESIZE):
    vf=np.zeros(SAMPLESIZE)
    tf=np.zeros(SAMPLESIZE)  
    ratio=int(round(TSAMPLING/tstep))
    ind=np.arange(0,int(np.floor(len(v)/ratio)))*ratio
    if len(ind)>SAMPLESIZE:
        ind=ind[0:SAMPLING]
    vf[0:len(ind)]=v[ind]
    tf[0:len(ind)]=t[ind]
    for k in range(len(ind),SAMPLESIZE):
        tf[k]=tf[k-1]+TSAMPLING
    return vf,tf

def Addnoise(vrms,v):
    for i in range(len(v)):
        v[i]=v[i]+gauss(0,vrms)
    return v

DISPLAY=0 # show plots if wanted


path=str(sys.argv[1]) # path to folder conating all the simulated traces and antenna positions
# array.dat: containing all antennas of an array = 10 000
# antpos.dat: containing just the antenna positions which were actually simulated after cone approximation
# -> find antenna position from array.dat in antpos.dat -> get the index ID from antpos.dat -> get out_ID.txt
# -> if antenna position from array.daz not fournd in antpos.dat -> create zero array with time trace starting at 0ns, 2ns,....

#units:
TSAMPLING=2e-9 #sec, eq to 500 MHz, sampling of system 
FREQMIN=50e6 #Hz  # frequencies which will be used in the later analysis: 50-200MHz
FREQMAX=200e6 #Hz, 250MHz
#tstep=t[1]-t[0]#1e-9, sec, time bins in simulations
#print " time binning sims ", tstep
SAMPLESIZE=int(3e-6/TSAMPLING) #=1500, 3e-6sec length # traces lenth of system
vrms=15 #uvolts, noise level



####Handing over one antenna or a whole array  ---- just handing over list of antennas possible at the moment
#if len(sys.argv)==3: # just one specif antenna handed over
    #start=int(sys.argv[2]) # antenna ID
    #end=start+1
    #print "single antenna with ID: ", str(start)," handed over"

if len(sys.argv)<3: # grep all antennas from the antenna file
  
    #positions=np.genfromtxt(path+'/antpos.dat')
    #start=0
    #end=len(positions)
    #print end
    
    positions=np.genfromtxt(path+'/array.dat') # full list of antennas
    pos=positions.tolist()
    start=0
    end=len(positions)
    print "number of antennas in array ", end
    
    positions_sim=np.genfromtxt(path+'/antpos.dat') # list of simulated antenns positions
    pos_sim=positions_sim.tolist()

    ###### loop  over l over all antennas in the full array list
    for l in range(start,end):
        print 'nr ', str(l)
        
## find the corresponding traces if antenna psoition was included in simulations        
        try:
            # compare array position with actual antenna positions, if position simulated read in the out_ID.txt file
            antenna_ID=pos_sim.index(pos[l]) # ID of simulated antenna 
            print "    found ", antenna_ID
            
            #NOTE: hand over path to file and antenna idea
            voltage_trace=path+"/out_"+str(antenna_ID)+".txt"
            print "try reading in ", voltage_trace
            try: 
                text=np.loadtxt(voltage_trace)#'out_128.txt')
            except IOError:
                print "IOError ..."
                continue
            t=text[:,0] # in s
            vx=text[:,1] #EW axis antenna
            vy=text[:,2] #NS axis antenna


            tstep=t[5]-t[4]#1e-9, sec, time bins in simulations
            #print "antenna ", str(l)," time binning sims ", tstep, 'total trace length :', t[-1]-t[0], len(t)

            
            if DISPLAY==1:
                plt.plot(t*1e9,vx) # s*1e9 = ns
                plt.plot(t*1e9,vy)
                plt.xlabel('Time [ns]')
                plt.ylabel('Voltage [uV]')
                plt.show()
                
            #Filtering in frequency band
            vx=Filtering(vx,tstep,FREQMIN,FREQMAX)
            vy=Filtering(vy,tstep,FREQMIN,FREQMAX)
            if DISPLAY==1:
                plt.plot(t*1e9,vx)
                plt.plot(t*1e9,vy)
                plt.xlabel('Time [ns]')
                plt.ylabel('Voltage [uV]')
                plt.show()
            
            
## fill up antenna which were not included in the sim            
        except ValueError:
            print "   antenna position not simulated" # empty time traces have to be created and just handed over to the rest
            
            tstep=1e-9# sec, time bins in simulations
            nbins=1024 # bin number from simulations
            
            
            vx=np.zeros(nbins)
            vy=np.zeros(nbins)
            t=np.fromfunction(lambda i: i*tstep, (nbins,), dtype=float)
            
            
            #Filtering in frequency band
            vx=Filtering(vx,tstep,FREQMIN,FREQMAX)
            vy=Filtering(vy,tstep,FREQMIN,FREQMAX)
  
            if DISPLAY==1:
                plt.plot(t*1e9,vx)
                plt.plot(t*1e9,vy)
                plt.xlabel('Time [ns]')
                plt.ylabel('Voltage [uV]')
                plt.show()
                
                
                
## create the fake traces

        #Filtering in frequency band
        vx=Filtering(vx,tstep,FREQMIN,FREQMAX)
        vy=Filtering(vy,tstep,FREQMIN,FREQMAX)
        if DISPLAY==1:
            plt.plot(t*1e9,vx)
            plt.plot(t*1e9,vy)
            plt.xlabel('Time [ns]')
            plt.ylabel('Voltage [uV]')
            plt.show()

        vx,tx=Digitization(vx,t,tstep,TSAMPLING,SAMPLESIZE)
        vy,ty=Digitization(vy,t,tstep,TSAMPLING,SAMPLESIZE)
        if DISPLAY==1:    
            plt.plot(vx)
            plt.plot(vy)
            plt.xlabel('Time [2 ns bins]')
            plt.ylabel('Voltage [uV]')
            plt.show()

        vx=Addnoise(vrms,vx)
        vy=Addnoise(vrms,vy)
        if DISPLAY==1:
            plt.plot(vx)
            plt.plot(vy)
            plt.xlabel('Time [2 ns bins]')
            plt.ylabel('Voltage [uV]')
            plt.show()


        outfile=path+"/fake_"+str(l)+".txt"
        f = file(outfile,"w")
        for i in np.arange(len(tx)):
                print >>f,"%e	%1.3e	%1.3e" % (tx[i], vx[i], vy[i] )  # time in s         
        f.close()
        print "saved as ", outfile
