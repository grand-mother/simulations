import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from random import *
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import sys
import time
import re

def Filtering(v,tstep,FREQMIN,FREQMAX):
    F=rfftfreq(len(v))/tstep #len(v) points between 0 and 1/2tstep
    V=rfft(v)
    V[F<FREQMIN]=0
    V[F>FREQMAX]=0
    return irfft(V)

def Digitization(v,t,tstep,TSAMPLING):
    ratio=int(round(TSAMPLING/tstep))
    ind = np.arange(0,len(v),ratio)
    vf=v[ind]
    tf=t[ind]
    return vf,tf

def Addnoise(vrms,v):
    return v+np.random.normal(0,vrms,size=np.shape(v))
    
def treat(path,id):
  #units:
  TSAMPLING=2e-9 #sec, eq to 500 MHz, sampling of system
  FREQMIN=50e6 #Hz  # frequencies which will be used in the later analysis: 50-200MHz
  FREQMAX=200e6 #Hz, 250MHz
  tstep=1e-9 #=t[1]-t[0], sec, time bins in Zhaires simulations
  #SAMPLESIZE=int(3e-6/TSAMPLING) #=1500 samples for 3musec length
  #SAMPLESIZE=1500 #=1500 samples for 3musec length
  SAMPLESIZE=300
  hw = SAMPLESIZE/2
  vrms=15 #uvolts

  #NOTE: hand over path to file and antenna idea
  voltage_trace=str(path)+"/out_"+str(id)+".txt"
  
  try:
    text=np.loadtxt(voltage_trace)
  except:
    return
  print voltage_trace,'now loaded.'
  t=text[:,0]
  t = t-np.min(t) # in s
  vx=text[:,1] #EW axis antenna
  vy=text[:,2] #NS axis antenna
  ## Filter
  vxf=Filtering(vx,tstep,FREQMIN,FREQMAX)
  vyf=Filtering(vy,tstep,FREQMIN,FREQMAX)
  txmax = np.argmin(vxf)
  vxf = vxf[np.argmin(vxf)-hw:np.argmin(vxf)+hw]
  txmax = np.argmin(vxf)
  vyf = vyf[np.argmin(vyf)-hw:np.argmin(vyf)+hw]
  tymax = np.argmin(vyf)
  t = t[txmax-hw:txmax+hw]
  vx = vx[txmax-hw:txmax+hw]
  vy = vy[tymax-hw:tymax+hw]
  
  tns = t*1e9
  print 'Raw signals pp amp =',max(vx)-min(vx),max(vy)-min(vy)
  print 'Filtered signals pp amp @ t=',txmax,tymax,':',max(vxf)-min(vxf),max(vyf)-min(vyf)
  if max(vxf)-min(vxf)<15 and max(vyf)-min(vyf)<15:
    return

  DISPLAY=0
  if DISPLAY==1:
      plt.figure(1)
      plt.plot(tns,vx) # s*1e9 = ns
      plt.plot(tns,vy)
      plt.xlabel('Time (ns)')
      plt.ylabel('Voltage ($\mu$V)')  
      plt.plot(tns,vxf,'b--')
      plt.plot(tns,vyf,'g--')
      plt.xlabel('Time (ns)')
      plt.ylabel('Filtered Voltage ($\mu$V)')
      plt.xlim(0,np.max(tns))
      imaxx = np.argmin(abs(vxf))
      imaxy = np.argmin(abs(vyf))
      #plt.show()

  # Add noise
  vxn=Addnoise(vrms,vxf*1e-6)
  vyn=Addnoise(vrms,vyf*1e-6)
  if DISPLAY==1:
      plt.figure(3)
      plt.subplot(211)
      plt.plot(tns,vxn)
      plt.ylabel('Voltage+Noise X ($\mu$V)')
      plt.xlim(0,np.max(tns))
      plt.subplot(212)
      plt.plot(tns,vyn)
      plt.xlim(0,np.max(tns))
      plt.xlabel('Time (ns)')
      plt.ylabel('Voltage+Noise Y ($\mu$V)')

  ## Digitize
  vxnd,tx=Digitization(vxn,t,tstep,TSAMPLING)
  vynd,ty=Digitization(vyn,t,tstep,TSAMPLING)
  vxnd = vxnd[vxnd!=0]
  vynd = vynd[vynd!=0]

  ## Oversample
  tck = interpolate.splrep(tx, vxnd, s=0)
  vxndo = interpolate.splev(t, tck, der=0)
  tck = interpolate.splrep(ty, vynd, s=0)
  vyndo = interpolate.splev(t, tck, der=0)
  vyndo[0:10]=0
  vxndo[0:10]=0
  vyndo[-10:]=0
  vxndo[-10:]=0

  if DISPLAY == 1:
    plt.figure(72)
    plt.subplot(211)
    plt.plot(tns,vxndo,'o-g',label='Oversd')
    plt.plot(tx*1e9,vxnd,'o-b',label='Digitized')
    plt.plot(tns,vxn,'r',label='Original')
    plt.plot((tns[txmax],tns[txmax]),(-max(abs(vxnd))*0.8,max(abs(vxnd))*1.2),'--r',lw=3)
    plt.ylabel('Signal ($\mu$V)')
    plt.subplot(212)
    plt.plot(tns,vyndo,'o-g',label='Oversd')
    plt.plot(tx*1e9,vynd,'o-b',label='Digitized')
    plt.plot(tns,vyn,'r',label='Original')
    plt.xlabel('Time (ns)')
    plt.ylabel('Signal ($\mu$V)')

  vxfd,tx=Digitization(vxf,t,tstep,TSAMPLING)
  vyfd,ty=Digitization(vyf,t,tstep,TSAMPLING)
  vxfd = vxfd[vxfd!=0]
  vyfd = vyfd[vyfd!=0]
  ppx = np.max(vxnd[txmax/2-10:txmax/2+10])-np.min(vxnd[txmax/2-10:txmax/2+10])
  ppy = np.max(vynd[tymax/2-10:tymax/2+10])-np.min(vynd[tymax/2-10:tymax/2+10])
  print 'Digitized signals pp amp:',ppx, ppy
  time.sleep(1)
  

  ## Now correlate to library
  #libpath = '/home/martineau/GRAND/GRAND/data/massProd/HS1orig/lib/'  #RM showers
  #libpath = '/home/martineau/GRAND/GRAND/data/CRs/showers/lib/'
  #libpath = '/home/martineau/GRAND/GRAND/data/CRs/showers/EE1E17.5-az0-zen100_evt0/split/'
  #libpath = path  #Analyse on itself
  #libpath = '/home/martineau/GRAND/GRAND/data/CRs/showers/EE1E19.5-az45-zen95_evt2/split/'
  #libpath = '/home/martineau/GRAND/GRAND/data/CRs/showers/EE1E17.5-az135-zen95_evt8/split'
  libpath = '/home/martineau/GRAND/GRAND/data/CRs/showers/EE1E18.0-az135-zen95_evtxx/split'
  xmaxc,xmaxnormc,xtmaxc = [],[],[]
  ymaxc,ymaxnormc,ytmaxc = [],[],[]
  sid = []
  for name in os.listdir(libpath):
    if not name.endswith("txt"):  # Not a voltage trace
	    continue
    filename = os.path.join(libpath, name)
    text=np.loadtxt(filename)#
    tl=text[:,0] # in s
    vxl=text[:,1] #EW axis antenna
    vyl=text[:,2] #NS axis antenna
    if np.argmin(vxl)<hw or np.argmin(vxl)>len(vxl)-hw:
      continue
    if np.argmin(vyl)<hw or np.argmin(vyl)>len(vyl)-hw:
      continue
    vxl = vxl[np.argmin(vxl)-hw:np.argmin(vxl)+hw]
    vyl = vyl[np.argmin(vyl)-hw:np.argmin(vyl)+hw]
    vxlf=Filtering(vxl,tstep,FREQMIN,FREQMAX)
    vylf=Filtering(vyl,tstep,FREQMIN,FREQMAX)
    stretchx = (np.max(vxndo)-np.min(vxndo))/(np.max(vxlf)-np.min(vxlf))
    stretchy = (np.max(vyndo)-np.min(vyndo))/(np.max(vylf)-np.min(vylf))   
    if (stretchx<0.1) or (stretchx>1000): # Don't stretch signals too much
      continue
    if (stretchy<0.1) or (stretchy>1000):
      continue
    vxlfdnorm = vxlf*stretchx
    vylfdnorm = vylf*stretchy
    corx = np.correlate(vxndo,vxlfdnorm,"full")
    cory = np.correlate(vyndo,vylfdnorm,"full")
    deltaoptx = len(vxlfdnorm)-np.argmin(vxlfdnorm)+txmax  # Expected position of correlation max when i=option = "full"
    xoffset = np.argmax(corx)-deltaoptx+1
    xmaxc.append(max(corx))
    xtmaxc.append(xoffset)
    nsigx = (max(corx)-np.mean(corx))/np.std(corx)
    xmaxnormc.append(nsigx)
    deltaopty = len(vylfdnorm)-np.argmin(vylfdnorm)+tymax  # Expected position of correlation max when i=option = "full"
    yoffset = np.argmax(cory)-deltaopty+1
    ymaxc.append(max(cory))
    ytmaxc.append(yoffset)
    nsigy = (max(cory)-np.mean(cory))/np.std(cory)
    ymaxnormc.append(nsigy)
    sid.append(name)
    #print deltaoptx,np.argmin(corx),xoffset,nsigx
    #print deltaopty,np.argmin(cory),yoffset,nsigy

 
    if DISPLAY==1:
      print 'Loading libray signal',filename
      plt.figure(41)
      plt.subplot(311)
      plt.plot(vxndo,'og-',label='Noisy')
      plt.xlim(0,len(vxndo))
      plt.plot(vxlfdnorm,'g:',lw=3,label='Original')
      plt.plot((txmax,txmax),(-max(abs(vxnd))*0.8,max(abs(vxnd))*1.2),'--r',lw=3)
      plt.grid(True)
      plt.subplot(312)
      plt.plot(corx,'k')
      plt.xlim(0,len(corx))
      plt.xlabel('Time (samples)')
      plt.grid(True)
      plt.subplot(313)
      plt.hist(corx,100)
 
      plt.figure(42)
      plt.subplot(311)
      plt.plot(vyndo,'og-',label='Noisy')
      plt.xlim(0,len(vyndo))
      plt.plot(vylfdnorm,'g:',lw=3,label='Original')
      plt.plot((tymax,tymax),(-max(abs(vynd))*0.8,max(abs(vynd))*1.2),'--r',lw=3)
      plt.grid(True)
      plt.subplot(312)
      plt.plot(cory,'k')
      plt.xlim(0,len(cory))
      plt.xlabel('Time (samples)')
      plt.grid(True)
      plt.subplot(313)
      plt.hist(cory,100)
     
      print deltaoptx,np.argmin(corx),xoffset,max(corx),nsigx
      print deltaopty,np.argmin(cory),yoffset,max(cory),nsigy

      plt.show()

  xmaxc = np.array(xmaxc)
  xtmaxc = np.array(xtmaxc)
  xmaxnormc = np.array(xmaxnormc)
  ymaxc = np.array(ymaxc)
  ytmaxc = np.array(ytmaxc)
  ymaxnormc = np.array(ymaxnormc)
  bestx = np.argmax(xmaxc)
  besty = np.argmax(ymaxc)
  print 'Best x:',xmaxc[bestx],xmaxnormc[bestx],xtmaxc[bestx],sid[bestx]
  print 'Best y:',ymaxc[besty],ymaxnormc[besty],ytmaxc[besty],sid[besty]

  plt.figure(12)
  plt.subplot(211)
  plt.plot(xtmaxc,xmaxc,'+')
  plt.title("X channel")
  plt.ylabel('Correlation Coef')
  plt.subplot(212)
  plt.plot(xtmaxc,xmaxnormc,'+')
  plt.xlabel('$\Delta$t truth (samples)')
  plt.ylabel('N StdDev')

  plt.figure(13)
  plt.subplot(211)
  plt.plot(ytmaxc,ymaxc,'+')
  plt.ylabel('Correlation Coef')
  plt.title("Y channel")
  plt.subplot(212)
  plt.plot(ytmaxc,ymaxnormc,'+')
  plt.xlabel('$\Delta$t truth (samples)')
  plt.ylabel('N StdDev')

  plt.figure(4)
  plt.subplot(211)
  plt.plot(tx*1e9,vxnd,'or-')
  plt.plot(tns,vxndo,'og-')
  # Load best
  filename = os.path.join(libpath, sid[bestx] )
  #print 'Loading',filename
  text=np.loadtxt(filename)
  vxl=text[:,1] #EW axis antenna
  vxl = vxl[np.argmin(vxl)-hw:np.argmin(vxl)+hw]
  vxlf=Filtering(vxl,tstep,FREQMIN,FREQMAX)
  vxlfnorm = vxlf*ppx/(np.max(vxlf)-np.min(vxlf))
  corx = np.correlate(vxndo,vxlfnorm,"full")
  shift = len(vxlfnorm)-np.argmax(corx)-1
  
  #deltaoptx = len(vxlfdnorm)-np.argmin(vxlfdnorm)+txmax  # Expected position of correlation max when i=option = "full"
  #xoffset = np.argmin(corx)-deltaoptx+1
  #print deltaoptx,np.argmin(corx),xoffset,shift    
  
  vxs = np.roll(vxlfnorm, -shift)
  if len(tns)==len(vxs):
    plt.plot(tns,vxs,'k',lw=2)
  plt.plot(tns,vxlfnorm,'k--')
  plt.ylabel('ChX signal ($\mu$V)')
  plt.subplot(212)
  plt.plot(ty*1e9,vynd,'og-')
  plt.plot(tns,vyndo,'or-')
  # Load best
  filename = os.path.join(libpath, sid[besty] )
  #print 'Loading',filename
  text=np.loadtxt(filename)
  vyl=text[:,2] #EW axis antenna
  vyl = vyl[np.argmin(vyl)-hw:np.argmin(vyl)+hw]
  vylf=Filtering(vyl,tstep,FREQMIN,FREQMAX)
  vylfnorm = vylf*ppy/(np.max(vylf)-np.min(vylf))
  cory = np.correlate(vyndo,vylfnorm,"full")
  shift = len(vylfnorm)-np.argmax(cory)-1
  vys = np.roll(vylfnorm, -shift)
  if len(tns)==len(vys):
    plt.plot(tns,vys,'k',lw=2)
  plt.plot(tns,vylfnorm,'k--')
  plt.xlabel('Time (ns)')
  plt.ylabel('ChY signal ($\mu$V)')
  
  #plt.show()
  
  bestIDx = int(re.split('\.|_',sid[bestx])[1])
  bestIDy = int(re.split('\.|_',sid[besty])[1])
  outfile="testCorrel.txt"
  f = file(outfile,"a")
  print >>f,"%3d %3.1f %3.1f %3.1f %3.1f %3d %3d %3.1f %3.1f %3d %3d %3.1f %3.1f" % (int(id),max(vxf)-min(vxf),max(vyf)-min(vyf),ppx,ppy,bestIDx,xtmaxc[bestx],xmaxc[bestx],xmaxnormc[bestx],bestIDy,ytmaxc[besty],ymaxc[besty],ymaxnormc[besty])
  f.close()

  #
  #f = file(outfile,"w")
  #for i in np.arange(len(tx)):
  #  print >>f,"%1.5e	  %1.2e   %1.2e" % (tx[i], vx[i], vy[i] )  # .2e as in the output of the sim
  #f.close()

if __name__ == '__main__':
   treat(sys.argv[1],sys.argv[2])
