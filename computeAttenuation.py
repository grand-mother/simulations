#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time
import modules
import numpy as np
DISPLAY = 0
if DISPLAY:
  import matplotlib.pyplot as pl
  import matplotlib.colors as colors
  import matplotlib.cm as cm
  from mpl_toolkits.mplot3d import Axes3D

RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
MAPDIR = "/home/martineau/GRAND/GRAND/data/maps/ASTER-GDEM2/"
sys.path.append(RETRODIR)
from grand_tour import Topography
from retro.event import EventIterator
#sys.path.append("/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/scripts/")

Rt = 6371000  # Earth radius (m)

def getGroundAltitudeFlat(x):
# Compute z componant of Earth surface at location (x[0],x[1]) in GRAND ref
   mod_cm = np.linalg.norm([x[0],x[1],x[2]+Rt])
   return Rt*((x[2]+Rt)/mod_cm-1)


def diffLossFlat(d,h,lam):
# Function to compute attenuation due to diffraction (implementation of ITU-R P.526-14) assuming Earth-curvature only along track
# Input parameters: 
# d: distance from source to antenna
# h: source height above ground
# lam: signal wavelength
# Output parameters:
# attenuation (dB)
# testx (validity test) 

    ae = 8500000 # Effective Earth radius
    hant = 5  # Antenna heigth
    dlos = np.sqrt(2*ae)*(np.sqrt(h)+np.sqrt(hant))  # Marginal distance along LoS
    #print 'd,dlos,h=',d,dlos,h
    if d<dlos:  # ==> use method 3.2
      #print "Distance =",d[i],"km<dlos =",dlos,", LoS not over the horizon. Computing attenuation."
      m = d*d/(4*ae*(h+hant))
      c = (h-hant)/(h+hant)
      b = 2*np.sqrt((m+1)/(3*m))*(np.cos(np.pi/3+np.arccos(3*c*np.sqrt((3*m)/pow(m+1,3))/2)/3))
      d1 = d*(1+b)/2
      d2 = d-d1
      h22 = ((h-d1*d1/(2*ae))*d2+(hant-d2*d2/(2*ae))*d1)/d  # Clearance height
      hreq = 0.552*np.sqrt(d1*d2*lam/d)
      #print "d1 =",d1,"m, d2 =",d2,"m"
      if h22>hreq:  # No attenuation
        #print "Clearence height =",h22,"m>hreq =",hreq,"m, no attenuation!"
        return 0, 0	
      else:  # Recompute effective Earth radius
	ae = 0.5*pow(d/(np.sqrt(h)+np.sqrt(hant)),2)
        #print "Clearence height =",h22/1e3,"<hreq =",hreq/1e3, ", associated effective Earth radius aem =",ae/1e3
 
    # Now compute attenuation
    y1 = 2*np.power(np.pi*np.pi/(lam*lam*ae),0.3333)*h
    if y1>2:
      g1 = 17.6*np.sqrt(y1-1.1)-5*np.log10(y1-1.1)-8
    else:  
      g1 = 20*np.log10(y1+0.1*np.power(y1,3))
    
    y2 = 2*np.power(np.pi*np.pi/(lam*lam*ae),0.3333)*hant
    g2 = 20*np.log10(y2+0.1*np.power(y2,3))

    x = np.power(np.pi/(lam*ae*ae),0.3333)*d
    if x<1.6:
        fx = -20*np.log10(x)-5.6488*np.power(x,1.425)
    else:
        fx = 11+10*np.log10(x)-17.6*x

    def delta(y):
      return 0.5*(1+np.tanh((0.5*np.log10(y)-0.255)/0.3))
    
    testx = x-np.sqrt(y1)*delta(y1)-np.sqrt(y2)*delta(y2)
    #if testx>1.096:
    #  print "**** Warning!!! Computation is not valid!!! testx =",testx
    
    dK = fx + g1 + g2   
    #print 'd,dlos=',d,dlos
    if d<dlos: 
      if h22<=hreq:
        dK = (1-h22/hreq)*dK
        return dK, testx
      else:
       print '*** Should not be here!!!'        	
    else:  # d>dlos
       return dK,testx
      

def compute_ray(r0, r1, lam): 
# Compute parameters for diffraction attenuation computation
# Input parameters
# r0: source location
# r1: antenna location
# lam: wavelength
# Output parameters:
# FresnelRange (m)
# z position @ Fresnel range location (x,y,z) in GRAND ref (m)
# Slope of plane traj between source and antenna (deg)
# Max offset to plane traj between source and antenna (m)
# RMS offset to plane traj between source and antenna (m)
    
    flat = False
    v = r1-r0
    if np.linalg.norm(v)>76000:  # Working only on antenna closer than 76km from Xmax
      print "Xmax too far"
      return None, None, None, None, None
    vn = v/np.linalg.norm(v)
    s = np.arange(0,np.linalg.norm(v),100)  # 100m step

    def compute_fresnelRadius(d1,lam):
    	d = max(d1)
    	d2 = d-d1
    	R = np.sqrt(d1*d2*lam/(d1+d2))
    	return R

    R = compute_fresnelRadius(s, lam)
    zt, zr, dalt = np.zeros(s.shape), np.zeros(s.shape), np.zeros(s.shape)
    for i, si in enumerate(s):         
        xi, yi, zi = r0 + si * vn	
        zr[i] = zi
	if flat == True:
	  zt[i]=getGroundAltitudeFlat([xi,yi,zi])
	  dalt[i] = np.linalg.norm([xi,yi,zi+Rt])-Rt
          #print i,xi,yi,zi,zt[i]
        else:
	  zt[i] = topo.ground_altitude(xi, yi)
          dalt[i] = zr[i]-zt[i]-20  # Ray height above ground (approx: ray height = CM-R_T in principle!!!)
    diffa = R>dalt # diffraction area

    if DISPLAY:
      pl.ion()
      pl.figure()
      pl.subplot(211)
      pl.plot(s, zt, "k")
      pl.plot(s, zr, "k--")
      pl.title('Antenna position {0} {1} {2}'.format(r1[0],r1[1],r1[2]))
      pl.subplot(212)
      pl.plot(s, -dalt, "k+-")
      pl.plot(s, -R, "b+-")
      pl.plot(s[diffa], -R[diffa], "r+-")
      pl.grid(True)
      pl.show() 
      raw_input()
	
    if np.sum(diffa)<1:  # Always flying
      #print 'Always flying'
      return 0, 0, -1000, 0, 0

    
    diffInd = np.where(R>dalt)[0]
    continuous =  np.sum(np.diff(diffInd)>3)<1 # Check if continuous - Allow up to one 3*100=300m discontinuity only
    #print diffInd,np.diff(diffInd)>2,continuous
    if continuous and (diffInd[-1]-len(R))<2:  # Continuous Fresnel range + ending at the antenna 
      fir = diffInd[0]  # First point of Fresnel range
      dFresnelRange = s[-1]-s[fir]  # Range where diffraction effects come into play: 
      aFlat = (zt[-1]-zt[fir])/dFresnelRange
      #print s[fir],zt[fir],s[las],zt[las],aFlat 
      y0 = zt[fir]-aFlat*s[fir]
      yFlat = aFlat*s[diffa]+y0
      dhmax = np.max(zt[diffa]-yFlat)  # Max deviation from plane ground
      dhrms = np.std(zt[diffa]-yFlat)  # Std deviation from plane ground
      sdeg = np.arctan2(zt[-1]-zt[fir],dFresnelRange)*180./np.pi 
      #print "slope=",sdeg
      return dFresnelRange, zt[fir], sdeg, dhmax, dhrms
    else:
      #print 'Distant or non continuous obstacle'
      return -1000, -1000, -1000, -1000, -1000  
 	 
	  
def fresnel(event):
    """Extract the relevant tau info from an event"""
    global origin, topo 
    flat = False  # Use topography of flat ground (with Earth curvature)
    DISPLAY = 0
    #50MHz: lam=6m
    #100MHz: lam=3m
    c0 = 299792458
    freq = np.arange(20,300)*1e6 # 100MHz
    lam = c0/freq
    
    if flat==False and event["origin"] != origin: #  Update the topography handle if required
      print "Loading new topo tile..."
      latitude, longitude = origin = event["origin"]
      print "(lat, lont)",latitude, longitude
      topo = Topography(latitude=latitude, longitude=longitude,path=MAPDIR, stack_size=121)     
      print 'Done.'
      
    # Get decay parameters
    _, _, r0, u, _, _ = event["tau_at_decay"]
    r0, u = map(np.array, (r0, u))
    tagall = event["tag"]
    tag = int(tagall.split(".")[-1])
    
    # Compute the shower energy
    shower_energy = 0.
    for (pid_, momentum) in event["decay"]:
    	aid = abs(pid_)
    	if aid in (12, 13, 14, 16):
    	    continue
    	shower_energy += sum(m**2 for m in momentum)**0.5
    zenith = np.arccos(-u[2])  # Zenith angle @ injection point (Zhaires convention, Radians)   

    # Get Xmax location
    Xmax_primary = modules._getXmax("pion", shower_energy,_) # Slant Xmax 
    #decayHeight = r0[2]-getGroundAltitudeFlat(r0)
    decayHeight = np.linalg.norm([r0[0],r0[1],r0[2]+Rt])-Rt
    Xmax_height, Xmax_distance = modules._dist_decay_Xmax(zenith, decayHeight, Xmax_primary) # d_prime: distance from decay point to Xmax
    r1 = r0 + u * Xmax_distance # Xmax position
    print 'Shower vector:',u
    print 'Injection point:',r0
    print 'Xmax:',r1
    print 'Xmax distance (km) =',Xmax_distance/1e3
    #raw_input()
    
    # Get antenna data
    ra = np.array(event["antennas"])[:, :3]
    
    attfile=tagall+".att"  
    if os.path.isfile(attfile):
      os.remove(attfile)   # Delete old file
    fid=open(attfile,'a')		

    ## Loop on antennas in the shower
    for i in range(len(ra[:, 1])):
      print '*** Antenna ',i,ra[i, :]

      dBAtt = np.ones((1,len(lam)))
      # Loop n all frequencies
      for j in range(len(lam)):
 	# Determine Fresnel range
 	r, f, a, m, s= compute_ray(r1,ra[i, :],lam[j]) # Ray is along line of sigth from Xmax to antenna
 	#print r, m, s, a
 
 	if r is None: # Antenna is to far
	  break  
	if r<0:  # Non continuous LoS ==> skip this frequency
	    continue
 	if r>0:
 	  # Compute attenuation in Fresnel range
 	  v = np.array(r1-ra[i,:])
 	  vn = v/np.linalg.norm(v)
 	  rx = ra[i,:]+vn*r  #Ray location at position of Fresnel range
 	  hx = rx[2]-f
 	  #print 'Source distance, height:',r,hx
 	  dbAtt2, ind2 = diffLossFlat(r,hx,lam[j])
 	if r==0:
 	  dbAtt2, ind2 = 0, 0

 	dBAtt[0,j] = dbAtt2
	
      dBAtt = np.insert(dBAtt,0,int(i))
      np.savetxt(fid, dBAtt[np.newaxis],fmt="%s",delimiter=' ',newline='\n')
      
def process(path):
    """Summarise a set of event files"""
    global origin, topo 
    origin, topo = None, None
    for name in os.listdir(path):
        if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        print "o Processing", filename
	t0 = time.time()
        for event in EventIterator(filename):
            fresnel(event)
        print "  --> Done in {:.1f} s".format(time.time() - t0)
    

if __name__ == "__main__":
    print "Usage: >python computeFresnel.py <path to json file>"
    process(sys.argv[1])
	
