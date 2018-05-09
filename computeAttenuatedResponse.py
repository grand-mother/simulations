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

DISPLAY = 0
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


def filt(f):
    #testFiltering()
    d = np.loadtxt(f)
    t = d[:,0]
    ns = d[:,1]
    ew = d[:,2]
    vert = d[:,3]
    t = (t-t[0])*1e9 #ns
    tstep=1e-9  #s
    fs = 1/tstep
    nsf = Filtering(ns,fs,FREQMIN,FREQMAX)
    ewf = Filtering(ew,fs,FREQMIN,FREQMAX)
    vertf = Filtering(vert,fs,FREQMIN,FREQMAX)
    
    return np.array([max(nsf)-min(nsf),max(ewf)-min(ewf),max(vertf)-min(vertf)])
    
    
    
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
       
       if DISPLAY: 
         col = ['b','g','r']   
	 t = (t-t[0])*1e9  
 	 pl.figure()
 	 pl.subplot(211)
 	 pl.plot(t,v,label='Free')
 	 pl.plot(t,vout[:,i-1],label='Attenuated')
 	 pl.xlabel("Time (ns)")
 	 pl.ylabel("Voltage ($\mu$V)")
 	 pl.xlim([0, max(t)])
 	 pl.subplot(212)
 	 pl.plot(f/1e6,np.abs(fv))
 	 pl.plot(f/1e6,np.abs(fva))
 	 pl.xlabel("Frequency (MHz)")
 	 pl.ylabel("FFT")
 	 pl.xlim([0, 300])
    
       # Now filter
       fs = 1/dt
       vout[:,i-1]  = Filtering(vout[:,i-1],fs,FREQMIN,FREQMAX)
       res = res+[t[imax],vout[imax,i-1],t[imin],vout[imin,i-1]]
    
       if DISPLAY:
	 print "Channel, Vpp",i-1,vout[imax,i-1]-vout[imin,i-1]
         pl.subplot(211)
         pl.plot(t,vout[:,i-1],label='Filtered')
	 pl.legend(loc='best')
	 pl.show()
	 raw_input()
	 
    vcom = vout[:,1]+vout[:,2] 
    imax = np.argmax(vcom,axis=0)
    imin = np.argmin(vcom,axis=0)
    res = res+[t[imax],vcom[imax],t[imin],vcom[imin]]
    return res


def process(jsonpath,attpath=None,tarpath=None):
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
    if attpath == None:
      attpath = jsonpath+"/"
    if tarpath == None:
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
 
    # Load attenuation
    print "Loading attenuation table",attf
    if os.path.isfile(attf):
      attt = np.loadtxt(attf)
      antid = attt[:,0]
      att = attt[:,1:]
    else:
      print "No file",attf,"! Abort."
      continue
 
    # Now loop on events in json file (Should only be one)
    for evt in EventIterator(jsonf):
      print "Now computing attenuated voltages"
      voltage=[]
      time_peaks=[]
      # Now loop on all antennas
      for i in range(int(max(antid))):
        
  	file_freespace = untardir+"/"+target+"/out_{0}.txt".format(i)
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
    shutil.rmtree(untardir)
   
if __name__ == "__main__":
    if len(sys.argv)==1:
    	print "Usage: >python computeAttenuatedResponse.py <path to json> [<path to att>] [<path to tgz>]"
    else:	
      print len(sys.argv)
      attpath, tarpath = None, None
      if len(sys.argv)>1:
        jsonpath = sys.argv[1]
      if len(sys.argv)>2:
        attpath = sys.argv[2]
      if len(sys.argv)>3:
        tarpath = sys.argv[3]
      
      process(jsonpath,attpath,tarpath)
