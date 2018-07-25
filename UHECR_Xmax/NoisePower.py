
#==============================================================================

# Script provided by Stijn Buitink adapted for GRAND Xmax reco.

#==============================================================================

import numpy as np

def Calc_Noise(zen_rot, az_rot, lowco, hico):
   
    coreasfile = 'fake_noise_voltage.txt'

    data=np.genfromtxt(coreasfile)
    dlength=data.shape[0]
    
    tstep = (data[1,0]-data[0,0]) # simulation sampling
    
    poldata=np.ndarray([dlength,2])

    XYZ=np.zeros([dlength,3])
    XYZ[:,0]=data[:,1]
    XYZ[:,1]=data[:,2]
    XYZ[:,2]=data[:,3]
    
    
    # Convert to on-sky coordinates (n, theta, phi) to prepare for application of antenna model
    poldata[:,0] = -1.0/np.sin(zen_rot)*XYZ[:,2]
    poldata[:,1] = np.sin(az_rot)*data[:,2] + np.cos(az_rot)*data[:,1]
    spec=np.fft.rfft(poldata, axis=-2)
    
    
    ### Substitute Antenna model by simple Dipol antenna with flat bandpass (50-350MHz) NOTE: add antenna response model here
    
    
    # Apply antenna model
    freqhi = 0.5/tstep/1e6 # MHz
    freqstep = freqhi/(dlength/2+1) # MHz   
        	
    instr_spec=np.ndarray([dlength/2+1,2],dtype=complex)    
    instr_spec[:,0] = spec[:,0]
    instr_spec[:,1] = spec[:,1]
    
    
    #Apply window (filtering to desired frequency band) and reduce maximum frequency to acquire downsampled signal
    fb = int(np.floor(lowco/freqstep))
    lb = int(np.floor(hico/freqstep)+1)

    ospow0=np.abs(spec[:,0])*np.abs(spec[:,0])  # filtered time series pol0
    ospow1=np.abs(spec[:,1])*np.abs(spec[:,1])  # filtered time series pol1
        
    filteredpower=np.sum(np.array([np.sum(ospow0[fb:lb+1]),np.sum(ospow1[fb:lb+1])])/(dlength/2.)*tstep)
    
    print(filteredpower)
    
    return filteredpower

