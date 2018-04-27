import sys
import glob
#sys.path.append("/home/martineau/GRAND/soft/neutrinos/simulations")
import modules
#import computeVoltage_HorAnt as comp
import os.path
import numpy as np
import pylab as pl
from scipy.signal import butter, lfilter
from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.interpolate import interp1d

sys.path.append("/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/scripts")
from retro.event import EventIterator, EventLogger
#sys.path.insert(0, '/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/scripts')
from plotShowerCaras import doScatterCol

from computeFresnel import fresnel

exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]
pl.style.use("/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/deps/mplstyle-l3/style/l3.mplstyle")

# Flat
#target ="/home/martineau/GRAND/GRAND/data/massProd/E.9e18_Z.89_A.334_La.41_Lo.87_H.663_D.42276570031824766"   
#target="/home/martineau/GRAND/GRAND/data/massProd/E.9e18_Z.90_A.283_La.42_Lo.86_H.1105_D.14135798811818680"
#target="/home/martineau/GRAND/GRAND/data/massProd/E.9e18_Z.92_A.115_La.42_Lo.87_H.2500_D.11451299271365943"
#target="/home/martineau/GRAND/GRAND/data/massProd/E.9e18_Z.95_A.10_La.41_Lo.86_H.4672_D.36379192182258783"  

path = "/home/martineau/GRAND/GRAND/data/massProd/HS1vertfix/"
# HS1
#target= "E.9e18_Z.92_A.157_La.42_Lo.87_H.2755_D.12509011140557780" 
#target= "E.9e18_Z.89_A.87_La.42_Lo.87_H.1648_D.29501657169930681" 
target= "E.9e18_Z.89_A.316_La.42_Lo.86_H.1497_D.6301798171336067"
#target= "E.7e17_Z.89_A.313_La.42_Lo.86_H.1501_D.21507668022063004"
jsonf = path+target+".voltage.json"
attf = path+target+".att"

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

def testFiltering():
    fs = 1e10  #10GHz
    frange = np.arange(1,30)*1e7  # Frequency range
    resp = np.zeros(frange.shape)
    t = np.arange(20000)*1e-10  # Time (s)
    for f in range(len(frange)):
      testx = np.sin(2*np.pi*frange[f]*t)
      testxf = Filtering(testx,fs,FREQMIN,FREQMAX)
      #pl.figure(f)
      #pl.plot(t*1e9,testx)
      #pl.plot(t*1e9,testxf)
      #pl.show()
      pos = np.arange(1000,len(t)) # Cut begining because unclean
      resp[f] = 20*np.log10((np.max(testxf[pos])-np.min(testxf[pos]))/(np.max(testx[pos])-np.min(testx[pos])))
      
    pl.figure()
    pl.plot(frange/1e6,resp)
    pl.xlabel("Frequency (MHz)")
    pl.ylabel("Filter response (dB)")
    pl.grid(True)
    pl.show()
    raw_input()


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
    
    if 0:
      pl.plot(t,nsf,'--b')
      pl.plot(t,ewf,'--g')
      pl.plot(t,vertf,'--r')
    
    return np.array([max(nsf)-min(nsf),max(ewf)-min(ewf),max(vertf)-min(vertf)])
    
    
def display(f):
    #s = f.split("_")
    #print s
    d = np.loadtxt(f)
    t = d[:,0]
    ns = d[:,1]
    ew = d[:,2]
    vert = d[:,3]
    t = (t-t[0])*1e9
    
    if 0:
      pl.figure()
      pl.plot(t,ns)
      pl.plot(t,ew)
      pl.plot(t,vert)
      pl.xlabel("Time (ns)")
      pl.ylabel("Voltage ($\mu$V)")
      pl.xlim([0, max(t)])
      pl.grid(True)

    return np.array([max(ns)-min(ns),max(ew)-min(ew),max(vert)-min(vert)])
    
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
        
    t = (t-t[0])*1e9
    vout = np.zeros(shape=(len(t),3))
    for i in range(1,4):
       v = d[:,i]
       fv = rfft(v)
       fva = fv*att
       vout[:,i-1] = irfft(fva)
       
       if 0:
 	 pl.figure()
 	 pl.subplot(211)
 	 pl.plot(t,v)
 	 pl.plot(t,vout[:,i])
 	 pl.xlabel("Time (ns)")
 	 pl.ylabel("Voltage ($\mu$V)")
 	 pl.xlim([0, max(t)])
 	 pl.subplot(212)
 	 pl.plot(f/1e6,np.abs(fv))
 	 pl.plot(f/1e6,np.abs(fva))
 	 pl.xlabel("Frequency (MHz)")
 	 pl.ylabel("FFT")
 	 pl.xlim([0, 300])
 	 raw_input()
     
    return np.array([max(vout[:,0])-min(vout[:,0]),max(vout[:,1])-min(vout[:,2]),max(vout[:,2])-min(vout[:,2])])

############################################################################################"
##### Main program
############################################################################################"

pl.ion()
DISPLAY = 0

# Load antennas
for evt in EventIterator(jsonf):  # Should only be one
    ants = np.array(evt["antennas"])
    xants = ants[:,0]
    yants = ants[:,1]
    zants = ants[:,2]
   
    #fresnel(evt)  # TBD if diffAtt.tmp not up to date

# Load attenuation
attf = np.loadtxt(attf)
antid = attf[:,0]
att = attf[:,1:]
res = []
for i in range(int(max(antid))):  # Allow 2000 antennas at most
  file_freespace = path+target+"/out_{0}_0.txt".format(i)  
  file_ground = path+target+"/out_{0}.txt".format(i)  
  if os.path.isfile(file_freespace) and os.path.isfile(file_ground):
    print '*** Antenna',i,xants[i],yants[i],zants[i]
    
    #if np.any(isout == i):  #Skip antenna for which no attenuation was computed
    #  print "Skip antenna",i
    #  continue
    #Vpp_0 = display(file_freespace)
    Vpp_0 = attenuate(file_freespace,att[i,:])
    
    #Vfpp_0 = filt(file_freespace)  # Skip filtering for now
    Vfpp_0 = Vpp_0
    
    g = np.loadtxt(file_ground)
    Vpp_g = display(file_ground)
    #Vfpp_g = filt(file_ground)
    Vfpp_g = Vpp_g

    #print Vfpp_0[0],Vfpp_0[1],Vfpp_0[2],Vfpp_g[0],Vfpp_g[1],Vfpp_g[2]
    #raw_input()
    res.append([xants[i],yants[i],zants[i],Vfpp_0[0],Vfpp_0[1],Vfpp_0[2],Vfpp_g[0],Vfpp_g[1],Vfpp_g[2]])
        
  else:
    print "No voltage for antenna",i  
 
 
res = np.array(res)
xall = res[:,0]
yall = res[:,1]
Af_ns = res[:,3]
Af_ew = res[:,4]
Af_v = res[:,5]
Ag_ns = res[:,6]
Ag_ew = res[:,7]
Ag_v = res[:,8]

doScatterCol(xall/1e3,yall/1e3,Af_ns,xlab='SN (km)',ylab='EW (km)',clab='AmpNS ($\mu V$)')
pl.title('Free space - NS')
doScatterCol(xall/1e3,yall/1e3,Af_ew,xlab='SN (km)',ylab='EW (km)',clab='AmpEW ($\mu V$)')
pl.title('Free space - EW')
doScatterCol(xall/1e3,yall/1e3,Af_v,xlab='SN (km)',ylab='EW (km)',clab='AmpVert ($\mu V$)')
pl.title('Free space - Vert')

doScatterCol(xall/1e3,yall/1e3,Ag_ns,xlab='SN (km)',ylab='EW (km)',clab='AmpNS ($\mu V$)')
pl.title('Ground - NS')
doScatterCol(xall/1e3,yall/1e3,Ag_ew,xlab='SN (km)',ylab='EW (km)',clab='AmpEW ($\mu V$)')
pl.title('Ground - EW')
doScatterCol(xall/1e3,yall/1e3,Ag_v,xlab='SN (km)',ylab='EW (km)',clab='AmpVert ($\mu V$)')
pl.title('Ground - Vert')

diffNS = (Af_ns-Ag_ns)/Ag_ns*100
diffEW = (Af_ew-Ag_ew)/Ag_ew*100
diffV = (Af_v-Ag_v)/Ag_v*100
    
doScatterCol(xall/1e3,yall/1e3,diffNS,figId=88,xlab='SN (km)',ylab='EW (km)',clab='DiffNS (\%)')
pl.title('Diff NS')
doScatterCol(xall/1e3,yall/1e3,diffEW,figId=89,xlab='SN (km)',ylab='EW (km)',clab='DiffEW (\%)')
pl.title('Diff EW')
doScatterCol(xall/1e3,yall/1e3,diffV,figId=90,xlab='SN (km)',ylab='EW (km)',clab='DiffVert (\%)')
pl.title('Diff Vert')

pl.figure()
pl.subplot(311)
pl.hist(diffNS,100)
print 'NS:',np.mean(diffNS),np.std(diffNS)
pl.subplot(312)
pl.hist(diffEW,100)
print 'EW:',np.mean(diffEW),np.std(diffEW)
pl.subplot(313)
pl.hist(diffV,100)
print 'Vert:',np.mean(diffV),np.std(diffV)

raw_input()
