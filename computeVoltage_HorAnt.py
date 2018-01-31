#!/usr/bin/env python
import os
from os.path import  join
import sys
import math
import numpy as np
#import pylab as pl

wkdir = './'

import linecache
from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.interpolate import interp1d
import retro
from retro.event import EventIterator, EventLogger
import modules

EARTH_RADIUS=6370949. #m
azstep=5 #step in azimuth in npy file
freqscale=1 #freq*2 if h/2 and sizeant/2
outputpower=0 #if wanted output is power
loaded=1 #if antenna is loaded or not in npy file
particle_list=[22.0, 11.0, -11.0, 111.0, 211.0, -211.0, 221.0] # 22:gamma, 11:e+-, 111:pi0, 211:pi+-, 211:eta

#impRLC R = 300;C = 6.5e-12;L = 1e-6; 20 300 MHz
RLp=np.array([0.536733768083299,   0.840010121593293,   1.200896057862110,   1.600090229038176,   2.006667705049151,   2.381373444983652,   2.685327039754095,2.890920915220645,   2.989008053352027,   2.988573924755322,   2.910298546933116,   2.778572562448605,   2.615623513015058,   2.438774302675680,2.260069404472798,   2.087095391770503,   1.924138291818527,   1.773244204213567,   1.635038639773870 ,  1.509302421351773,   1.395352352338537,1.292281470704730,   1.199104504595046 ,  1.114841938895556,   1.038565549999183,   0.969420412709374 ,  0.906632958509609 ,  0.849511074453546, 0.797439919119660,   0.749875669027360,   0.706338496130072 ,  0.666405514438127,   0.629704091546769,   0.595905715891539 ,  0.564720490529909, 0.535892256306097,   0.509194310914832,   0.484425672932856 ,  0.461407833493318,   0.439981938133439,   0.420006344534532 ,  0.401354506652400, 0.383913141085314,   0.367580636870446 ,  0.352265674931463 ,  0.337886027975219 ,  0.324367515703664,   0.311643093771139 ,  0.299652058008190, 0.288339348095085 ,  0.277654937150177,   0.267553295648130,   0.257992919745857,   0.248935915510479 ,  0.240347631749694,   0.232196335171802, 0.22445292247744])
RLp=RLp*100
XLp=np.array([1.149833973427903 ,  1.347001618559050 ,  1.469876468210024,   1.496656923296413,   1.411838459831799,   1.213746617081856 ,  0.919238711558534, 0.561550538777911 ,  0.181259529550333 , -0.184790883266832 , -0.510938360781749,  -0.784380139061162 , -1.002688474655982 , -1.169915727151229, -1.293172262455534 , -1.380332931202411 , -1.438786540600538,  -1.474902574702375,  -1.493909155794964 , -1.499971154708314,  -1.496345170687206,  -1.485548051254960 , -1.469510769217091 , -1.449707994034062,  -1.427262501557595 , -1.403027192021063 , -1.377648559710691,  -1.351615388245273,  -1.325295941574338 , -1.298966315222550 , -1.272832045980507,  -1.247044599699961 , -1.221713972961578 , -1.196918345352935 , -1.172711490165158, -1.149128477825457 , -1.126190075642858 , -1.103906149182129 , -1.082278296775352,  -1.061301893203182,  -1.040967676805738,  -1.021262982755671,  -1.002172701363368 , -0.983680022166382,  -0.965767010753354 , -0.948415054722766 , -0.931605207084645 , -0.915318449184856 , -0.899535890421292,  -0.884238918293781 , -0.869409309431790 , -0.855029309984293 , -0.841081691988703 , -0.827549790949401 , -0.814417528766047 , -0.801669425292115, -0.789290601124629])
XLp=XLp*100
fr=np.arange(20,301,5)

# Load antenna response files
fileleff_x=wkdir+'HorizonAntenna_SNarm_leff_loaded.npy' # 'HorizonAntenna_leff_notloaded.npy' if loaded=0, EW component, used for NS
fileleff_y=wkdir+'HorizonAntenna_EWarm_leff_loaded.npy' # 'HorizonAntenna_leff_notloaded.npy' if loaded=0, EW component, used for EW
fileleff_z=wkdir+'HorizonAntenna_Zarm_leff_loaded.npy' # 'HorizonAntenna_leff_notloaded.npy' if loaded=0, EW component, used for vert
freq1,realimp1,reactance1,theta1,phi1,lefftheta1,leffphi1,phasetheta1,phasephi1=np.load(fileleff_x) ### this line cost 6-7s
RL1=interp1d(fr, RLp, bounds_error=False, fill_value=0.0)(freq1[:,0])
XL1=interp1d(fr, XLp, bounds_error=False, fill_value=0.0)(freq1[:,0])
freq2,realimp2,reactance2,theta2,phi2,lefftheta2,leffphi2,phasetheta2,phasephi2=np.load(fileleff_y) ### this line cost 6-7s
RL2=interp1d(fr, RLp, bounds_error=False, fill_value=0.0)(freq2[:,0])
XL2=interp1d(fr, XLp, bounds_error=False, fill_value=0.0)(freq2[:,0])
freq3,realimp3,reactance3,theta3,phi3,lefftheta3,leffphi3,phasetheta3,phasephi3=np.load(fileleff_z) ### this line cost 6-7s
RL3=interp1d(fr, RLp, bounds_error=False, fill_value=0.0)(freq3[:,0])
XL3=interp1d(fr, XLp, bounds_error=False, fill_value=0.0)(freq3[:,0])

# #===========================================================================================================
# def GRANDtoNEC(zenith=None, azimuth =None):
# #===========================================================================================================
#     zen = (180-zenith)
#     azim = azimuth + 90 # az_ant=az_GRAND +90deg
#     if azim>360:
#       azim = azim-360
#
#     return [zen,azim]
#
# #===========================================================================================================
# def NECtoGRAND(zenith=None, azimuth =None):
# #===========================================================================================================
#     zen = (180-zenith)
#     azim = azimuth - 90 # az_GRAND=az_ant -90deg
#     if azim>360:
#         azim = azim-360
#     elif azim<0:
#         azim = azim+360
#
#     return [zen,azim]

#===========================================================================================================
def get_voltage(time1, Ex, Ey, Ez, ush=[1, 0, 0], alpha=0, beta=0, typ="X"):
#===========================================================================================================
    # Note: azim & zenith are in GRAND convention

    def TopoToAntenna(u,alpha,beta): #from coordinates in the topography frame to coordinates in the antenna
        alpha=alpha*np.pi/180 #around y
        beta=beta*np.pi/180 #around x
        cb = np.cos(beta)
        sb = np.sin(beta)
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        # rotx = np.array([[1,0,0],[0,cb,-sb],[0,sb,cb]])
        # rotx = np.linalg.inv(rotx)  # Referential rotates ==> use inverse matrix
        roty = np.array([[ca,0,sa],[0,1,0],[-sa,0,ca]])
        #roty = np.linalg.inv(roty)  # Not to be used since y is facing backwards in (x,y,z)
        rotz = np.array([[cb,-sb,0],[sb,cb,0],[0,0,1]])
        rotz = np.linalg.inv(rotz)
        rotyz=roty.dot(rotz)  # beta and then alpha rotation. This induces a EW component for x arm

        # Now rotate along zp so that we are back with x along NS
        xarm = [1,0,0]  #Xarm
        xarmp = rotyz.dot(xarm)  # Rotate Xarm along slope
        # Compute antrot, angle vs NS, and then rotate back along NS (angle = -antrot)
        antrot = math.atan2(xarmp[1],xarmp[0])*180/np.pi
        #print "antrot = ",antrot
        cz = np.cos(antrot*np.pi/180)
        sz = np.sin(antrot*np.pi/180)
        rotzant = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
        rotzant = np.linalg.inv(rotzant)
        rottot = rotzant.dot(rotyz)

        [xp,yp,zp] = rottot.dot(u)
        return np.array([xp,yp,zp])

    # Load proper antenna response matrix
    if typ=="X":
       fileleff = fileleff_x
       freq=freq1
       realimp=realimp1
       reactance=reactance1
       theta=theta1
       phi=phi1
       lefftheta=lefftheta1
       leffphi=leffphi1
       phasetheta=phasetheta1
       phasephi=phasephi1
       RL=RL1
       XL=XL1
    if typ=="Y":
       fileleff = fileleff_y
       freq=freq2
       realimp=realimp2
       reactance=reactance2
       theta=theta2
       phi=phi2
       lefftheta=lefftheta2
       leffphi=leffphi2
       phasetheta=phasetheta2
       phasephi=phasephi2
       RL=RL2
       XL=XL2
    if typ=="Z":
       fileleff = fileleff_z
       freq=freq3
       realimp=realimp3
       reactance=reactance3
       theta=theta3
       phi=phi3
       lefftheta=lefftheta3
       leffphi=leffphi3
       phasetheta=phasetheta3
       phasephi=phasephi3
       RL=RL3
       XL=XL3

    # Compute effective theta, phi in antenna tilted frame (taking slope into account, with x=SN)
    ushp = TopoToAntenna(ush,alpha,beta)  # Xmax vector in antenna frame
    zen=np.arccos(ushp[2])*180/np.pi  # Zenith in antenna frame
    azim=math.atan2(ushp[1],ushp[0])*180/np.pi
    if azim>360:
        azim = azim-360
    elif azim<0:
        azim = azim+360
    # print ush, ushp, alpha, beta
    # print zen,azim
    # print [1,0,0],TopoToAntenna([1,0,0],alpha,beta)
    # print [0,1,0],TopoToAntenna([0,1,0],alpha,beta)
    # print [0,0,1],TopoToAntenna([0,0,1],alpha,beta)
    # alpha()

    if typ=='X':
        print "Zenith & azimuth in antenna framework:",zen, azim
    if zen>90:
        print "Signal originates below antenna horizon! No antenna response computed. Abort."
        return([],[])
    # Now take care of Efield signals
    delt = time1[1]-time1[0];
    Fs = 1/delt
    timeoff=time1[0] # time offset, to get absolute time
    time1 = (time1-time1[0]) #reset to zero
    # Rotate Efield to antenna frame (x along actual arm)
    Etot=np.array([Ex,Ey,Ez])
    [Exp,Eyp,Ezp] = TopoToAntenna(Etot,alpha,beta)
    szen = np.sin(zen*np.pi/180);
    czen = np.cos(zen*np.pi/180);
    saz = np.sin(azim*np.pi/180);
    caz = np.cos(azim*np.pi/180);
    #amplituder = szen*(caz*Exp+saz*Eyp)+czen*Ezp
    amplitudet = czen*(caz*Exp+saz*Eyp)-szen*Ezp
    amplitudep = -saz*Exp+caz*Eyp
    # if typ == "Z":
    #     pl.figure(12)
    #     pl.plot(Exp)
    #     pl.plot(Eyp)
    #     pl.plot(Ezp)

    ##################################
    ### all the settings for the 3 different antenna arms:

    nfreq=len(freq[:,0])
    f=np.zeros(nfreq)
    RA=np.zeros(nfreq)
    XA=np.zeros(nfreq)
    ltr=np.zeros(nfreq)
    lta=np.zeros(nfreq)
    lpr=np.zeros(nfreq)
    lpa=np.zeros(nfreq)
    if azstep==5:
        roundazimuth=round(azim/10)*10+round((azim-10*round(azim/10))/5)*5
    elif azstep==1:
        roundazimuth=round(azim)
    else:
        print('Error on azimuth step!')
        return(0)
    if roundazimuth>=91 and roundazimuth<=180:
        roundazimuth=180-roundazimuth
    if roundazimuth>=181 and roundazimuth<=270:
        roundazimuth=roundazimuth-180
    if roundazimuth>=271 and roundazimuth<=360:
        roundazimuth=360-roundazimuth

    for i in range(nfreq):
        f[i]=freq[i,0]*freqscale
        indtheta=np.nonzero(theta[i,:]==round(zen))[0]
        indphi=np.nonzero(phi[i,:]==roundazimuth)[0]
        indcom=np.intersect1d(indtheta,indphi)
        ltr[i]=lefftheta[i,indcom]
       	lta[i]=np.deg2rad(phasetheta[i,indcom]) #*np.pi/180
       	lpr[i]=leffphi[i,indcom]
       	lpa[i]=np.deg2rad(phasephi[i,indcom]) #*np.pi/180
        if loaded==0:
            RA[i]=realimp[i,0]
            XA[i]=reactance[i,0]
            Rlefft=ltr[i]*np.cos(lta[i])
            Xlefft=ltr[i]*np.sin(lta[i])
            Rleffp=lpr[i]*np.cos(lpa[i])
            Xleffp=lpr[i]*np.sin(lpa[i])
            Rleqt=((Rlefft*RL[i]-Xlefft*XL[i])*(RA[i]+RL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            Xleqt=((Rlefft*RL[i]+Xlefft*XL[i])*(XA[i]+XL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            ltr[i]=np.sqrt(Rleqt**2+Xleqt**2)
            lta[i]=np.arccos(Rleqt/ltr[i])
            Rleqp=((Rleffp*RL[i]-Xleffp*XL[i])*(RA[i]+RL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            Xleqp=((Rleffp*RL[i]+Xleffp*XL[i])*(XA[i]+XL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            lpr[i]=np.sqrt(Rleqp**2+Xleqp**2)
            print(Rleqp,lpr[i])
            lpa[i]=np.arccos(Rleqp/lpr[i])

    if loaded==0:#phases are not unwrap! so:
        for i in range(1,nfreq):
            while lpa[i]-lpa[i-1]<-np.pi: #180*np.pi/180:
                lpa[i]=lpa[i]+ 2.*np.pi #360*np.pi/180
            while lpa[i]-lpa[i-1]> np.pi: #180*np.pi/180:
                lpa[i]=lpa[i]- 2.*np.pi #360*np.pi/180
            while lta[i]-lta[i-1]<-np.pi: #180*np.pi/180:
                lta[i]=lta[i]+ 2.*np.pi #360*np.pi/180
            while lta[i]-lta[i-1]>np.pi: #180*np.pi/180:
                lta[i]=lta[i]- 2.*np.pi #360*np.pi/180

    #print(round(zenith),roundazimuth,f,ltr,lta,lpr,lpa,Fs)
    ###############################

    fmin=f[0]
    fmax=f[-1]
    f=f*1e6

    nf  = int(2**np.floor(np.log(len(amplitudet))/np.log(2)))
    while Fs/nf > fmin*1e6:   # <== Make sure that the DFT resolution is at least fmin.
        nf *= 2
    F = rfftfreq(nf)*Fs

    modulust = interp1d(f, ltr, bounds_error=False, fill_value=0.0)(F)
    phaset   = interp1d(f, lta, bounds_error=False, fill_value=0.0)(F)
    modulusp = interp1d(f, lpr, bounds_error=False, fill_value=0.0)(F)
    phasep   = interp1d(f, lpa, bounds_error=False, fill_value=0.0)(F)
    if outputpower:
        RLinter  = interp1d(f, RL, bounds_error=False, fill_value=0.0)(F)

    phaset -= phaset[0] # Switch the phase origin to be consistent with a real signal.
    phasep -= phasep[0] # Switch the phase origin to be consistent with a real signal.

    #if we want P=V2/RL -> incorrect
    if outputpower:
        modulusp[RLinter!=0]=modulusp[RLinter!=0]/np.sqrt(RLinter[RLinter!=0])
        modulust[RLinter!=0]=modulust[RLinter!=0]/np.sqrt(RLinter[RLinter!=0])

    #B and D are V in freq domain, they are complex
    A = rfft(amplitudet, nf)
    ct = np.cos(phaset)
    st = np.sin(phaset)
    B = np.zeros(A.shape)
    B[1:-1:2] = modulust[1:-1:2]*(A[1:-1:2]*ct[1:-1:2]-A[2:-1:2]*st[2:-1:2])
    B[2:-1:2] = modulust[2:-1:2]*(A[1:-1:2]*st[1:-1:2]+A[2:-1:2]*ct[2:-1:2])
    B[0]  = A[0]*modulust[0]
    B[-1] = A[-1]*modulust[-1]

    C = rfft(amplitudep, nf)
    cp = np.cos(phasep)
    sp = np.sin(phasep)
    D = np.zeros(C.shape)
    D[1:-1:2] = modulusp[1:-1:2]*(C[1:-1:2]*cp[1:-1:2]-C[2:-1:2]*sp[2:-1:2])
    D[2:-1:2] = modulusp[2:-1:2]*(C[1:-1:2]*sp[1:-1:2]+C[2:-1:2]*cp[2:-1:2])
    D[0]  = C[0]*modulusp[0]
    D[-1] = C[-1]*modulusp[-1]

    #we should apply 1/sqrt(RL) to the real part and then put vt and vp squared, after ifft
    vt=irfft(B)
    vp=irfft(D)

    if outputpower:
        vt=vt**2
        vp=vp**2

    voltage = vp + vt
    timet     = np.arange(0, len(vt))/Fs
    timep     = np.arange(0, len(vp))/Fs

    #print '    Peak to peak voltage amplitude = ', max(voltage) - min(voltage),'muV'

    return(voltage, timet+timeoff)


#===========================================================================================================
def inputfromjson(path,json_file):
#===========================================================================================================
    # shower you are interested in
    showerID = str(path.split('/')[-1])
    if not showerID:
        showerID = str(path.split('/')[-2])

    # find that shower in the json file
    event = [evt for evt in EventIterator(json_file) if evt["tag"]==showerID][0]

    ### DECAY
    decay_pos=event["tau_at_decay"][1]
    #print "decay position: ", decay_pos
    injection_height=decay_pos[2]
    decay_pos=decay_pos+np.array([0.,0.,EARTH_RADIUS]) # corrected for earth radius
    #print "decay position after correction: ", decay_pos
    decay_altitude=event["tau_at_decay"][3]
    #print "decay decay_altitude: ", decay_altitude

    ### ANGLES
    v=event["tau_at_decay"][2]# shower direction, assuming decay products strongly forward beamed
    zenith_sim = np.degrees(np.arccos(np.dot(v, decay_pos) / np.linalg.norm(decay_pos))) # zenith in GRAND conv.
    #print "theta: ", zenith_sim
    #orthogonal projection of v onto flat plane to get the azimuth
    x=np.array([1.,0.,0.]) #NS
    y=np.array([0.,1.,0.]) #EW
    proj_v= np.dot(v,x)*x + np.dot(v,y)*y
    azimuth_sim = np.degrees(np.arccos(np.dot(proj_v, x))) # azimuth in GRAND conv., rt NORTH
    if proj_v[1]<0.: # y component of projection negativ, means azimuth >180deg
        azimuth_sim = 360.-azimuth_sim
    #print "azimuth: ", azimuth_sim

    ### ENERGY
    ep_array=np.zeros(len(event["decay"])-1)
    for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
        if float(event["decay"][i][0]) in particle_list: # just for valid particles
            pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1]
            ep_array[i-1]=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
    #    print "particle ", str(i), "PID:",event["decay"][i][0]," energy in EeV: ", ep_array[i-1]*1e-9 #np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)* 1.e-9
    energy= np.sum(ep_array)* 1.e-9 # GeV in EeV
    #print "energy in EeV: ", energy

    ### PID primary
    #part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
    particle=int(np.argmax(ep_array) +1) # not forget about the inital neutrino and array start with 0
    PID= float(event["decay"][int(np.argmax(ep_array) +1)][0])
    el_list=[22.0, 11.0, -11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11':positron,  '111.0': 'pi0'
    if PID in el_list:
        primary="electron"
    else: # pion-like
        primary="pion"

    return zenith_sim,azimuth_sim,energy,injection_height,primary

#===========================================================================================================
def inputfromtxt(input_file_path):
#===========================================================================================================
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

    energy = energy *1e-18

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

#===========================================================================================================
def compute(opt_input,path, effective,zenith_sim, azimuth_sim, energy, injection_height, primary,json_file=None):
#===========================================================================================================

    if opt_input=='json':
        # shower you are interested in
        showerID = str(path.split('/')[-1])
        if not showerID:
            showerID = str(path.split('/')[-2])
        # Find that shower in the json file
        event = [evt for evt in EventIterator(json_file) if evt["tag"]==showerID][0]
        log_event = EventLogger(path=json_file)

    voltage=[]
    time_peaks=[]
    print "Zenith, azimuth=",zenith_sim, azimuth_sim

    ##########################################################################################
    ###Handing over one antenna or a whole array
    if opt_input=='txt':
        if len(sys.argv)>=6: # just one specif antenna handed over
            start=int(sys.argv[5]) # antenna ID
            end=start+1
        #    print "single antenna with ID: ", str(start)," handed over"
        if  len(sys.argv)<6: # grep all antennas from the antenna file
            positions=np.genfromtxt(path+'/antpos.dat')
            start=0
            end=len(positions)
        #    print "Array with ", end, " antennas handed over"
    elif opt_input=='json':
        if len(sys.argv)>=6: # just one specif antenna handed over
            start=int(sys.argv[5]) # antenna ID
            end=start+1
        #    print "single antenna with ID: ", str(start)," handed over"
        if  len(sys.argv)<6: # grep all antennas from the antenna file
            positions=np.array(event["antennas"],dtype=float)
            decay_pos=event["tau_at_decay"][1]
            positions = positions - [decay_pos[0],decay_pos[1],0.]
            #positions=np.genfromtxt(path+'/antpos.dat')
            start=0
            end=len(positions)
            #print "Array with ", end, " antennas handed over"
    elif opt_input=='manual':
        if len(sys.argv)>=10: # just one specif antenna handed over
            start=int(sys.argv[9]) # antenna ID
            end=start+1
        #    print "single antenna with ID: ", str(start)," handed over"
        if  len(sys.argv)<10: # grep all antennas from the antenna file
            positions=np.genfromtxt(path+'/antpos.dat')
            start=0
            end=len(positions)
        #    print "Array with ", end, " antennas handed over"

    if effective==1: # effective zenith caclculation needs Xmax position as input
    #    print "effective zenith calculated - Xmax position approximated ..."
        # Then compute Xmax
        caz = np.cos(np.deg2rad(azimuth_sim))
        saz = np.sin(np.deg2rad(azimuth_sim))
        czen = np.cos(np.deg2rad(zenith_sim))
        szen = np.sin(np.deg2rad(zenith_sim))

        Xmax_primary = modules._getXmax(primary, energy, np.deg2rad(zenith_sim)) # approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
        Xmax_height, Xmax_distance = modules._dist_decay_Xmax(np.deg2rad(zenith_sim), injection_height, Xmax_primary) # d_prime: distance from decay point to Xmax
        Xmax = Xmax_distance*np.array([caz*szen, saz*szen, czen])+np.array([0,0,injection_height])
    #   print 'Xmax=',Xmax_primary,' Xmax height=',Xmax_height,' Xmax distance =',Xmax_distance,'Xmax position= ',Xmax
    #   print 'Now computing Xmax position from injection height=',injection_height,'m and (zen,azim) values.'

    ###### loop  over l --- LOOP OVER ANTENNA ARRAY
    for l in range(start,end):
        efieldtxt=path+'/a'+str(l)+'.trace'
    #    print 'Wave direction: zenith = ', zenith_sim, ' deg, azimuth = ', azimuth_sim, 'deg. (GRAND conventions), mountain slope: ', alpha_sim, 'deg.'
    #    print 'Efield file: ', efieldtxt

        # Model the input signal.
        try:
            time1_sim, Ex_sim, Ey_sim,Ez_sim = np.loadtxt(efieldtxt,delimiter=' ',usecols=(0,1,2,3),unpack=True)
        except IOError:
            continue

        # NOTE: adapt to your time from whatever to s
        time1_sim= time1_sim*1e-9 # time has to be handed in s

        #print 'Now computing antenna response...'
        if effective==1:

            # Compute effective zenith
            # First get antenna position

            #print 'Reading antenna position from parameter input.'
            if (opt_input=='json' or opt_input=='txt') and (len(sys.argv)==11) :
                    x_sim = float(sys.argv[6])
                    y_sim = float(sys.argv[7])
                    z_sim = float(sys.argv[8])

                    # include a mountain slope - correction of zenith angle
                    alpha_sim=float(sys.argv[9])
                    beta_sim=float(sys.argv[10])

            elif (opt_input=='manual') and (len(sys.argv)==15) :
                    x_sim = float(sys.argv[10])
                    y_sim = float(sys.argv[11])
                    z_sim = float(sys.argv[12])

                    # include a mountain slope - correction of zenith angle
                    alpha_sim=float(sys.argv[13])
                    beta_sim=float(sys.argv[14])
            else :
                try :
                    if opt_input=='json':
                        x_sim,y_sim,z_sim = positions[l] #,alpha_sim,beta_sim
                        alpha_sim = 0.
                        beta_sim = 0.
                    else:
                        #print 'Trying to read antenna position from antpos.dat file...'
                        numberline = int(l) + 1
                        line = linecache.getline(path+'/antpos.dat', numberline)
                        #[x_sim, y_sim, z_sim] = map(float, line.split())
                        [x_sim, y_sim, z_sim, alpha_sim, beta_sim] = map(float, line.split())
                    #print 'Read antenna position from antpos.dat file... Antenna',l,' at position [', x_sim, y_sim, z_sim,'].'
                except :
                    print 'No antenna position file found, please put antpos.dat in', path, 'or enter check antenna informations in json file or enter antenna positions as arguments.'
                    sys.exit()

        Xant = [x_sim, y_sim, z_sim]
 	# Hack OMH 24/01
	#alpha_sim=10
	#beta_sim=0
        #Xant = [400000, 0 , 0]
        ush = Xmax-Xant
        ush = ush/np.linalg.norm(ush)  # Unitary vector pointing to Xmax from antenna pos
        voltage_NS, timeNS  = get_voltage( time1=time1_sim,Ex=Ex_sim, Ey=Ey_sim, Ez=Ez_sim, ush=ush, alpha=alpha_sim, beta=beta_sim, typ="X")
        voltage_EW, timeEW  = get_voltage( time1=time1_sim,Ex=Ex_sim, Ey=Ey_sim, Ez=Ez_sim, ush=ush, alpha=alpha_sim, beta=beta_sim, typ="Y")
        voltage_vert, timevert  = get_voltage( time1=time1_sim,Ex=Ex_sim, Ey=Ey_sim, Ez=Ez_sim, ush=ush, alpha=alpha_sim, beta=beta_sim, typ="Z")

        #pl.savetxt(path+'out_'+str(l)+'.txt', (timeEW, voltage_EW, voltage_NS), newline='\r\n')#, voltage_NS)) # is not working correctly
        if np.size(timeEW)>0:   # Dat was computed
          f = file(path+'/out_'+str(l)+'.txt',"w")
          print "OUTFILE : ", path+'/out_'+str(l)+'.txt'
          for i in np.arange(len(timeEW)):
            print >>f,"%1.5e	%1.2e	%1.2e	%1.2e" % (timeNS[i], voltage_NS[i], voltage_EW[i], voltage_vert[i] ) # same number of digits as input
          f.close()

        ###plots
        DISPLAY=0
        if DISPLAY==1:
            import pylab as pl
            import matplotlib.pyplot as plt
            plt.figure(1,  facecolor='w', edgecolor='k')
            plt.subplot(211)
            plt.plot(time1_sim*1e9,Ey_sim, label="Ey = EW")
            plt.plot(time1_sim*1e9,Ex_sim, label="Ex = NS")
            plt.plot(time1_sim*1e9,Ez_sim, label="Ez = UP")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Electric field (muV/m)')
            plt.legend(loc='best')
            plt.subplot(212)
            plt.plot(timeEW*1e9,voltage_EW, label="EW")
            plt.plot(timeEW*1e9,voltage_NS, label="NS")
            plt.plot(timeEW*1e9,voltage_vert, label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage (muV)')
            plt.legend(loc='best')

            plt.show()

##################################################################
##################################################################
        if opt_input=='json':
            #### additional output needed for later study, added in the json file
            # p2p voltage:  antenna ID, p2p EW, NS, UP, EW+NS
            voltage_com=np.copy(voltage_EW)
            for i in range (0, len(voltage_EW)):
                                voltage_com[i]+=voltage_NS[i]
            v_list =( str(l),  max(voltage_EW) - min(voltage_EW), max(voltage_NS) - min(voltage_NS), max(voltage_vert) - min(voltage_vert), max(voltage_com) - min(voltage_com)   )
            voltage.append( v_list )

            # time of peaks and value: t_EW_max, v_EW_max, t_EW_min, v_EW_min,.... EW, NS, vert, EW+NS
            import operator
            EW_ind_max, value = max(enumerate(voltage_EW), key=operator.itemgetter(1))
            EW_ind_min, value = min(enumerate(voltage_EW), key=operator.itemgetter(1))

            NS_ind_max, value = max(enumerate(voltage_NS), key=operator.itemgetter(1))
            NS_ind_min, value = min(enumerate(voltage_NS), key=operator.itemgetter(1))

            vert_ind_max, value = max(enumerate(voltage_vert), key=operator.itemgetter(1))
            vert_ind_min, value = min(enumerate(voltage_vert), key=operator.itemgetter(1))

            com_ind_max, value = max(enumerate(voltage_com), key=operator.itemgetter(1))
            com_ind_min, value = min(enumerate(voltage_com), key=operator.itemgetter(1))

            time_peaks.append( (round(timeEW[EW_ind_max],11),  voltage_EW[EW_ind_max], round(timeEW[EW_ind_min],11), voltage_EW[EW_ind_min],
                                round(timeNS[NS_ind_max],11), voltage_NS[NS_ind_max], round(timeNS[NS_ind_min],11), voltage_NS[NS_ind_min],
                                round(timevert[vert_ind_max],11), voltage_vert[vert_ind_max], round(timevert[vert_ind_min],11), voltage_vert[vert_ind_min],
                                round(timeEW[com_ind_max],11), voltage_com[com_ind_max], round(timeEW[com_ind_min],11), voltage_com[com_ind_min] )  )

############### end of loop over antennas
    if opt_input=='json':
        if len(voltage)==0:
            print "- effective zenith not fulfilled - NO VOLTAGE COMPUTED"
            log_event(**event)
        else:
            # add the additional informations to the shower event
            event['voltage'] = voltage # muV
            event['time_peaks'] = time_peaks # s, muV
            log_event(**event)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#===========================================================================================================
# Compute the time dependent voltage
#===========================================================================================================
if __name__ == '__main__':

    if len(sys.argv)<5:
        print """\
        Wrong minimum number of arguments. All angles are to be expressed in degrees and in GRAND convention.
        Usage:
        if json file or ZHAireS inp file input:
            python computevoltage.py [input option] [effective 0/1] [path to traces] [json file/inp file] [opt: antenna ID] [opt: antenna x,y,z,alpha,beta]
            example: python computevoltage.py json/txt 0/1 ./ ../Danton/*.json 7 100 100 1000 10 5

        if manual input:
            python computevoltage.py [input option] [effective 0/1] [path to traces] [primary] [zenith] [azimuth] [energy in EeV] [injection height above sea level in m] [opt: antenna ID] [opt: antenna x,y,z,alpha,beta]
            example: python computeVoltage_massProd.py manual 0/1 ./ proton 85 205 0.5 1500 7 100 100 1000 10 5
            example: python computeVoltage_massProd.py txt 1 ./split/ ./ZhairesShower.inp

        """
        ## -> computes voltage traces for EW, NS and Vertical antenna component and saves the voltage traces in out_'.txt (same folder as a'.trace)
        ## -> produces a new json file with copying the original one, but saves as well additional informations as p2p-voltages, and peak times and values in *.voltage.json in the same folder as the original json file
        sys.exit(0)

    print "READING INPUT PARAMETERS"

    # Decide where to retrieve the shower parameters : #json for json file, txt for ZHAireS input file or manual to hand them over by hand
    opt_input = str(sys.argv[1])
    print opt_input

    # decide if the effectice zenith should be calculated (1) or not (0)
    effective = float(sys.argv[2])
    print effective

    # which efield trace do you wanna read in. to be consistent the script works with the antenna ID
    path=sys.argv[3] #folder containing the traces and where the output should go to

    if opt_input=='txt':
        # Read the ZHAireS input (.inp) file to extract the primary type, the energy, the injection height and the direction
        inp_file = str(sys.argv[4])
        zenith_sim,azimuth_sim,energy,injection_height,primary = inputfromtxt(inp_file)
        json_file = None

    elif opt_input=='json':
        # Read the json file to extract the primary type, the energy, the injection height, and the direction
        json_file = str(sys.argv[4])
        zenith_sim,azimuth_sim,energy,injection_height,primary = inputfromjson(path,json_file)

    elif opt_input=='manual':
        primary = str(sys.argv[4])
        zenith_sim = float(sys.argv[5]) #deg
        azimuth_sim = float(sys.argv[6]) #deg
        energy = float(sys.argv[7]) #EeV
        injection_height = float(sys.argv[8]) #m above sea level
        json_file = None

    #print 'shower = ',zenith_sim,azimuth_sim,energy
    print "VOLTAGE COMPUTATION STARTED"

    compute(opt_input,path, effective, zenith_sim, azimuth_sim, energy, injection_height, primary,json_file)

    print "VOLTAGE COMPUTED"
