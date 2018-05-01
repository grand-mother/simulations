import sys
import glob
import time
import modules
import os.path
import tarfile
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from scipy.fftpack import rfft, irfft, rfftfreq            

CC = 1
if CC==1:
  RETRODIR = "/pbs/throng/trend/soft/sim/GRANDsim/retro/"
else:
  RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
sys.path.append(RETRODIR)
sys.path.append(RETRODIR+"lib/python/")

sys.path.append("/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/scripts")
from retro.event import EventIterator, EventLogger

exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]



def attenuate(f,attdB):
    
    # Compute attenuation coefs
    fint = np.arange(20,300)*1e6
    attDB = np.array(attdB)
    attLin = pow(10,attdB/20.)
    attLin[attdB==+1] = 0  # attdB == 1 corresponds to values for which attenuation could not be computed (usually, no direct LoS)
    
    # Now load data and compute fft
    d = np.loadtxt(f)
    t = d[:,0]
    dt = np.mean(np.diff(t))  # Compute time step
    fmax = 1/dt/2 
    f = np.linspace(0,fmax,len(t))  # Frequency range
    # Interpolate attenuation coefs at desired frequency values
    fatt = interp1d(fint,attLin)
    att = np.zeros(np.shape(f))
    sel = np.where( (f>30e6) & (f<299e6) )
    att[sel] = fatt(f[sel])
        
    #t = (t-t[0])*1e9
    vout = np.zeros(shape=(len(t),3))
    res = []
    for i in range(1,4):
       v = d[:,i]
       fv = rfft(v)
       fva = fv*att
       vout[:,i-1] = irfft(fva)
       
       imax = np.argmax(vout[:,i-1],axis=0)
       imin = np.argmin(vout[:,i-1],axis=0)
       res = res+[t[imax],vout[imax,i-1],t[imin],vout[imin,i-1]]
       
       if 0:      
 	 pl.figure()
 	 pl.subplot(211)
 	 pl.plot(t,v)
 	 pl.plot(t,vout[:,i-1])
 	 pl.xlabel("Time (ns)")
 	 pl.ylabel("Voltage ($\mu$V)")
 	 pl.xlim([0, max(t)])
 	 pl.subplot(212)
 	 pl.plot(f/1e6,np.abs(fv))
 	 pl.plot(f/1e6,np.abs(fva))
 	 pl.xlabel("Frequency (MHz)")
 	 pl.ylabel("FFT")
 	 pl.xlim([0, 300])
	 
    vcom = vout[:,1]+vout[:,2] 
    imax = np.argmax(vcom,axis=0)
    imin = np.argmin(vcom,axis=0)
    res = res+[t[imax],vcom[imax],t[imin],vcom[imin]]
    return res


def process(target):

  print "o Processing", target
  t0 = time.time()

  # Set target path
  #evtpath = "/home/martineau/GRAND/GRAND/data/massProd/flat/"
  evtpath = './'
  jsonf = evtpath+target+".voltage.json"
  attf = evtpath+target+".att"
  tarf =  evtpath+target+".tgz"

  if os.path.isfile(tarf):
    # Unpack tar file with traces
    print "Extracting",tarf
    tar = tarfile.open(tarf, "r:gz")
    tar.extractall(path=evtpath)
    tar.close()
    print "Done."
  else:
    print "No file",tarf,"! Abort."
    return
    
  # Load attenuation
  print "Loading attenuation table",attf
  if os.path.isfile(attf):
    attt = np.loadtxt(attf)
    antid = attt[:,0]
    att = attt[:,1:]
  else:
    print "No file",attf,"! Abort."
    return
  
  # Now loop on events in json file (Should only be one)
  for evt in EventIterator(jsonf):  
    print "Now computing attenuated voltages"
    voltage=[]
    time_peaks=[]
    # Now loop on all antennas
    for i in range(int(max(antid))):
 
      file_freespace = evtpath+target+"/out_{0}.txt".format(i)
      if os.path.isfile(file_freespace):  # Voltage was computed
    	res = attenuate(file_freespace,att[i,:])
    	v_list = (i,round(res[1]-res[3],3),round(res[5]-res[7],3),round(res[9]-res[11],3),round(res[13]-res[15],3))
    	voltage.append( v_list )
    	time_peaks.append(res)
 
    # Now write output to file
    print "Now logging results to",jsonf
    log_event = EventLogger(path=jsonf)
    # Add the additional informations to the shower event
    evt['voltage'] = voltage # muV
    evt['time_peaks'] = time_peaks # s, muV
    log_event(**evt)

  print "  --> Done in {:.1f} s".format(time.time() - t0)


   
if __name__ == "__main__":
    print "Usage: >python computeAttenuatedResponse.py <ID of target event>"
    process(sys.argv[1])
