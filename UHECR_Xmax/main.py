
#==============================================================================
# UHECR reconstruction for GRAND
#==============================================================================

import sys
import numpy as np
import matplotlib.pyplot as plt

import process_GRAND_Zhaires as process
import analyze_GRAND_Zhaires as analyze


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

#==============================================================================
# Script for filtering etc.
#==============================================================================

process.Process_Data(datadir, fileno, outputfile, files, lowco, hico, filename)

#==============================================================================
# Script for analysis
#==============================================================================

for i in np.arange(int(files)-1):
    print "\nSimulated event: %i" %i
    results = analyze.reverseAnalysis(outputfile, eventno=fileno, eventno2=str(txt[i]), simevent=i, outputfolder=outputf, SKAmap=Map)
    print "Final result: reconstr, real Xmax "
    print results [0:2]    
    r2[i]=np.array([results[0], results[1], results[2], np.abs(results[0]- results[1]), np.abs(results[0]- results[2]), results[3], results[7]])


#file = open("ECR190_0_500m/Histtest_CR190_0.dat", "r")
#r2 = np.loadtxt(file)

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


lim=200.
Xreco_new = r2.T[0][ (r2.T[0]-r2.T[1] < lim)&(r2.T[0]-r2.T[1] > -lim)]
Xreal_new = r2.T[1][(r2.T[0]-r2.T[1] < lim)&(r2.T[0]-r2.T[1] > -lim)]

mean_tot=np.mean(r2.T[0]-r2.T[1])
std_tot=np.std(r2.T[0]-r2.T[1])

mean_new=np.mean(Xreco_new-Xreal_new)
std_new=np.std(Xreco_new-Xreal_new)

max_diff = np.max(r2.T[3])
int_max_diff = int(max_diff/100+1)*100


print("\n------FINAL RESULTS------")

print("\nMean: %.3e" %mean_tot)
print("Std: %.3e" %std_tot)

print("\nMean new: %.3e" %mean_new)
print("Std new: %.3e" %std_new)
print("Typical 68%: %.3e" %position2)

print("\nmax(Xreco-Xreal): %e \n" %max_diff)


#==============================================================================
# Plot the reconstruction uncertainty as a histogram
#==============================================================================

fig=plt.figure(1,figsize=(8,6)) 

bins=np.arange(-int_max_diff-10,int_max_diff+10,10)
plt.hist(r2.T[0]-r2.T[1], bins=bins, histtype='step')
plt.hist(Xreco_new-Xreal_new, bins=bins, histtype='step')

#bins=np.arange(0,100,5)
#plt.hist(r2.T[3], bins=bins, histtype='step')
#plt.axvline(position2,color='k', linestyle='--')

plt.ylabel("Nr. of Sim.")
#plt.xlabel("|Xreco-Xreal| (g/cm$^2$)")
plt.xlabel("Xreco-Xreal (g/cm$^2$)")
plt.xlim(-int_max_diff-10,int_max_diff+10)

plt.show()
name = outputf+'Histtest2_{0}_unc{1}.pdf'.format(sys.argv[2], std_new)
plt.savefig(name)


#==============================================================================
# Save all of it as a text file for later analysis
#==============================================================================

namefile = outputf+'Histtest_{0}.dat'.format(sys.argv[2])
file= open(namefile, 'w')  
#file.write('rec, real Xmax, best sim xmax, rec-real, rec-best, xmaxreco_stijn, primary\n')
for s in range(0,int(files)):
	  file.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(r2[s,0],r2[s,1],r2[s,2],r2[s,3],r2[s,4],r2[s,5],r2[s,6]))
file.close()

