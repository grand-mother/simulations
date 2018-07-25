import numpy as np


def noisepower(bandwidthmin, bandwidthmax): # power in RCL circuit
    dnu=10
    
    folder="./noisepower/" # power in RCL circuit

    noisefile=folder+'GRAND_EW'+str(bandwidthmin)+'.txt'
    text=np.loadtxt(noisefile)
    lst=text[:,0]
    powertot=text[:,1] #V^2

    for i in range(bandwidthmin+dnu,bandwidthmax,dnu):
        noisefile=folder+'GRAND_EW'+str(i)+'.txt'
        text=np.loadtxt(noisefile)
        powertot=powertot+text[:,1]
    
    return powertot


def noisevoltage(bandwidthmin, bandwidthmax):#V^2
    dnu=10

    folder="./noisevoltage/" #V^2

    noisefile=folder+'GRAND_EW'+str(bandwidthmin)+'.txt'
    text=np.loadtxt(noisefile)
    lst=text[:,0]
    powertot=text[:,1] #V^2

    for i in range(bandwidthmin+dnu,bandwidthmax,dnu):
        noisefile=folder+'GRAND_EW'+str(i)+'.txt'
        text=np.loadtxt(noisefile)
        powertot=powertot+text[:,1]
    
    return powertot

freqmin=50
freqmax=200

voltagetot=noisevoltage(freqmin, freqmax)
print "in ", freqmin, " - ", freqmax, "MHz band :::  max V2: ", max(voltagetot), " mean V2: ", np.mean(voltagetot), " Vrms: ",np.sqrt(np.mean(voltagetot)*2.)
        
powertot=noisepower(freqmin, freqmax)
print "in ", freqmin, " - ", freqmax, "MHz band :::  max power: ", max(powertot), " mean power: ", np.mean(powertot)