import sys
from sys import argv

import numpy as np
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
import ownprocess_wolora4_Zhaires as ow
import process_func as prf
import ska_stijn_wolora4_Zhaires as ska
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Parabola model
def model(x, a, b, c):
    return a*(x-b)**2 +c

script = sys.argv[0]  ### name of script
datadir = sys.argv[1]  ### path to data read in
fileno= sys.argv[2] ### run nr
outputfile= sys.argv[3] ### Outputfile, Pickle

# Frequencies: 30-80 MHz or 50-200 MHz
lowco= int(sys.argv[5])
hico=int(sys.argv[6])

# Output file, pickle
outputf=sys.argv[7]

# Antenna file GRAND
Map = sys.argv[8]

# Handing over the simulation which should be used
filename = sys.argv[9]

l_sim=open(filename).read().splitlines()
txt = np.array(l_sim)
files = len(txt)

r=np.zeros([files,7])
r2=np.zeros([files,7])
r3=np.zeros([files,7])

print "\nOutput file: " + outputfile
print "datadir: " + datadir
print "fileno: " + fileno
print "outputf: " + outputf
print "Map: " + Map
print "filename: " + filename


## Filtering 50-200 MHz etc and writing pickle file;
## BUT no antenna model (used to project x,y,z component to reality, needed e.g. for triggering, polarisation study)
## TimeBin 1.25ns

# Script for filtering etc.
ow.Process_Data(datadir, fileno, outputfile, files, lowco, hico, filename)


for i in np.arange(int(files)-1):
    print "\nSimulated event: %i" %i
    results = ska.reverseAnalysis(outputfile, eventno=fileno, eventno2=str(txt[i]), simevent=i, outputfolder=outputf, SKAmap=Map)
    print "Final result: reconstr, real Xmax "
    print results [0:2]    
    r2[i]=np.array([results[0], results[1], results[2], np.abs(results[0]- results[1]), np.abs(results[0]- results[2]), results[3], results[7]])


a2 = r2[:,3]
a2.sort()

total=files

proz2=0
position2=0
for m in np.arange(int(files)):  # sort them
    if(proz2 < 0.68*total): # get the 68% reconstruction uncertainty
      proz2+=1
    if(proz2 >= 0.68*total):
      position2 = a2[m]
      break

m=np.mean(r2.T[3])
std=np.std(r2.T[3])
print("\n------FINAL RESULTS------") 
print("Std: %e" %std)
print("3-sigma: %e" %(3.*std))
print("5-sigma: %e" %(5.*std))

 
#### Plot the reconstruction uncertainty as a histogram: new fit
fig=plt.figure(1,figsize=(8,6)) 
 
plt.hist(r2.T[3], bins=20, histtype='step')
plt.axvline(position2,color='k', linestyle='--')

plt.ylabel("Nr. of Sim.")
plt.xlabel("|Xreco-Xreal| (g/cm$^2$)")

plt.show()
name = outputf+'Histtest2_{0}_unc{1}.pdf'.format(sys.argv[2], position2)
plt.savefig(name)

### Save all of it as a text file for later analysis
namefile = outputf+'Histtest_{0}.dat'.format(sys.argv[2])
file= open(namefile, 'w')  
#file.write('rec, real Xmax, best sim xmax, rec-real, rec-best, xmaxreco_stijn, primary\n')
for s in range(0,int(files)):
	  file.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(r2[s,0],r2[s,1],r2[s,2],r2[s,3],r2[s,4],r2[s,5],r2[s,6]))
file.close()

