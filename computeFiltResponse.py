import sys
import glob
import shutil
import time
import modules
import os.path
import tarfile
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter
from scipy.fftpack import rfft, irfft, rfftfreq            

DISPLAY=False
CC = 0
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


FREQMIN = 50e6
FREQMAX = 200e6

def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    y = lfilter(b, a, data)
    return y

def Filtering(x,fs,lowcut,highcut):
    return butter_bandpass_filter(x, lowcut, highcut, fs, order=5)


def filt(f=None,d=None):
# Filters signals (passed as argument or from file)
# Returns Vpp

    if d is None:  # No array of data
      if f is None:  # No file is given
        print "No data input! Abort."
	return None
      else: # Load file
        d = np.loadtxt(f)
     
    # Now filter
    t = d[:,0]
    dt = np.mean(np.diff(t))  # Compute time step
    fs = 1/dt
    v = np.array(d[:,1:])  # Raw signal
    nCh = np.shape(v)[1]
    vout = np.zeros(shape=(len(t),nCh))	
    res = []
    for i in range(nCh):
      vi = v[:,i]
      vout[:,i]  = Filtering(vi,fs,FREQMIN,FREQMAX)
      imax = np.argmax(vout[:,i],axis=0)
      imin = np.argmin(vout[:,i],axis=0)
      res = res+[t[imax],vout[imax,i],t[imin],vout[imin,i]]
 
      if DISPLAY:
        print "Channel, Vpp filtered:",i,vout[imax,i]-vout[imin,i]
        pl.subplot(211)
        tdis = (t-t[0])*1e9 
	pl.plot(tdis,vout[:,i],label='Filtered')
        pl.legend(loc='best')
        pl.show()
        raw_input()
    
    if nCh==3:  # We have X,Y,Z channels == we can compute X+Y response
      # Combined X+Y channel
      vcom = vout[:,1]+vout[:,2] 
      imax = np.argmax(vcom,axis=0)
      imin = np.argmin(vcom,axis=0)
      res = res+[t[imax],vcom[imax],t[imin],vcom[imin]]
      
    return np.array(res)

    
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
    v = np.array(d[:,1:])  # Raw signal
        
    vatt = np.zeros(shape=(len(t),3))
    res = []
    for i in range(3):
       vi = v[:,i]
       fv = rfft(vi)
       fva = fv*att
       vatt[:,i] = irfft(fva)
       
       if DISPLAY: 
         imax = np.argmax(vatt[:,i],axis=0)
         imin = np.argmin(vatt[:,i],axis=0)
	 print "Channel, Vpp attenuated:",i,vatt[imax,i]-vatt[imin,i]
         col = ['b','g','r']   
	 tdis = (t-t[0])*1e9  
 	 pl.figure()
 	 pl.subplot(211)
 	 pl.plot(tdis,vi,label='Free')
 	 pl.plot(tdis,vatt[:,i],label='Attenuated')
 	 pl.xlabel("Time (ns)")
 	 pl.ylabel("Voltage ($\mu$V)")
 	 pl.xlim([0, max(tdis)])
 	 pl.subplot(212)
 	 pl.plot(f/1e6,np.abs(fv))
 	 pl.plot(f/1e6,np.abs(fva))
 	 pl.xlabel("Frequency (MHz)")
 	 pl.ylabel("FFT")
 	 pl.xlim([0, 300])
       
       # Now filter signal
       dout = np.array([t]+[vatt[:,i]])  
       dout = np.transpose(dout)     
       resi = filt(d=dout)
       res = np.concatenate((np.array(res),np.array(resi)),axis=0)
    
    # Combined X+Y channel
    vcom = vatt[:,1]+vatt[:,2] 
    dout = np.array([t]+[vcom])
    dout = np.transpose(dout)   
    resi = filt(d=dout)
    res = np.concatenate((res,resi),axis=0)
    
    return res


def process(jsonpath,attpath=None,tarpath=None):
  doAtt = False
  untardir = "./tmp"
  try:
    os.stat(untardir)
  except:
    os.mkdir(untardir)   
  
  for name in os.listdir(jsonpath):
    if not name.endswith("json"):
    	continue
    jsonf = name
    target = jsonf[:-13]
  
    print "o Processing", target
    t0 = time.time()
    if attpath is None:
      attpath = jsonpath+"/"
    if tarpath is None:
      tarpath = jsonpath+"/"
    jsonf = jsonpath+"/"+jsonf 
    attf = attpath+target+".att"
    tarf =  tarpath+target+".tgz"

    if os.path.isfile(tarf):
      # Unpack tar file with traces
      print "Extracting",tarf
      tar = tarfile.open(tarf, "r:gz")
      tar.extractall(path=untardir)
      tar.close()
      print "Done."
    else:
      print "No file",tarf,"! Abort."
      continue
 
    if doAtt:
      # Load attenuation
      print "Loading attenuation table",attf
      if os.path.isfile(attf):
        attt = np.loadtxt(attf)
        antid = attt[:,0]
        att = attt[:,1:]
        print "Done."
      else:
        print "No file",attf,"! Abort."
        #continue
    else:
        print "No attenuation considered in this treatment because doAtt=",doAtt
	
    # Now loop on events in json file (Should only be one)
    for evt in EventIterator(jsonf):
      ants = np.array(evt["antennas"])
      nAnts = len(ants)
      
      print "Now computing attenuated + filtered voltages"
      voltage=[]
      time_peaks=[]
      # Now loop on all antennas
      for i in range(nAnts):
  	file_freespace = untardir+"/"+target+"/out_{0}.txt".format(i)
	if os.path.isfile(file_freespace):  # Voltage file exists
	  if doAtt:
	    res = attenuate(file_freespace,att[i,:])
	  else:
	    res = filt(file_freespace)
	  
	  res[np.isnan(res)] = 0;
	  v_list = (i,round(res[1]-res[3],3),round(res[5]-res[7],3),round(res[9]-res[11],3),round(res[13]-res[15],3))
	  voltage.append( v_list )
  	  time_peaks.append(list(res))
 
      # Now write output to file
      print "Now logging results to",jsonf
      log_event = EventLogger(path=jsonf)
      # Add the additional informations to the shower event
      evt['voltage'] = voltage # muV
      evt['time_peaks'] = time_peaks # s, muV
      log_event(**evt)

    print "  --> Done in {:.1f} s".format(time.time() - t0)
    shutil.rmtree(untardir)
   
if __name__ == "__main__":
    if len(sys.argv)==1:
    	print "Usage: >python computeAttenuatedResponse.py <path to json> [<path to att>] [<path to tgz>]"
    else:	
      attpath, tarpath = None, None
      if len(sys.argv)>1:
        jsonpath = sys.argv[1]
      if len(sys.argv)>2:
        attpath = sys.argv[2]
      if len(sys.argv)>3:
        tarpath = sys.argv[3]
      
      process(jsonpath,attpath,tarpath)
