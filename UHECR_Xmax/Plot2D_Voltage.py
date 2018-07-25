

import sys
from sys import argv
import numpy as np
from numpy import *

import matplotlib.pyplot as plt
#import plotdata 
import pylab
import os
from matplotlib.pyplot import cm 

from scipy.signal import hilbert
import operator

from mpl_toolkits.mplot3d import Axes3D

def mag(x):
    return numpy.sqrt(x.dot(x))


####################################
decay=[67407.46, 0.00, 6332.68] 


# python wavefront.py pathtosim/
path= sys.argv[1]


## NOTE: as far as I remember you were able to reconstruct the correct antenna positions for antpos.dat-file using Nicolas script. they should be the sam one a sgiven in the simualtion 
## inpute file
posfile = path +'antpos.dat' # 
positions=np.genfromtxt(posfile)# positions[:,0]:along North-South, positions[:,1]: along East-West ,positions[:,2]: Up, in m    
x_pos= positions.T[0]
y_pos= positions.T[1]
z_pos= positions.T[2]
#x_pos, y_pos, z_pos= np.loadtxt(posfile,delimiter=' ',usecols=(0,1,2),unpack=True)# positions[:,0]:along North-South, positions[:,1]: along East-West ,positions[:,2]: Up, in m 

number_ant= len(x_pos) # number of positions you have
print(number_ant)


peaktime= np.zeros(number_ant)
peakamp= np.zeros(number_ant)
peakamp2= np.zeros(number_ant)
for i in range(0, number_ant): # loop over all antennas
    try: 
        #txt=np.loadtxt(path+ 'a'+str(i)+'.trace') # just pure simulation, in the ideal case you read in the voltage trace you get after applying computevoltage.py
    ## txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all in muV/m
    
        txt=np.loadtxt(path+ 'out_'+str(i)+'.txt')
    
    
    # now it depends of which parameter we need the peak amplitude: the pure trace, the hilbert envelope, combind signal of NS...
    # - eg. hilbert(txt.T[1]) would calculate the hilbert envelope of the NS component
    # - peak mpltitde of one component eg: abs(txt.T[1])
    # - combining NS and EW amplitude = np.sqrt(txt.T[1]**2. + txt.T[2]**2.)
        amplitude= np.sqrt(txt.T[1]**2.  + txt.T[2]**2.) ## x,y or NS, EW
    
        index, peak_amp = max(enumerate(amplitude), key=operator.itemgetter(1)) # get the index and the value of the peakamp
    #print(txt.T[0,index], peak_amp)
        peaktime[i]=txt.T[0,index] # get the time of the peak amplitude
        peakamp[i]=amplitude[index] # get the time of the peak amplitude
    except:
        peaktime[i]= 0.
        peakamp[i]= 0.
        

           
       
#diff_amp=peakamp-peakamp2    

#print peaktime    

## NOTE: for other experiments you know where the shower hits the ground (shower core) from the information from the particle detectors. In GRAND we do not have that. The questions is 
## how can we get it :) Is it like a fitting of function and look for the minimum.... This is what Frank hopefully tells us or do you already have an idea ...

x_pos=x_pos+decay[0]
y_pos=y_pos+decay[1]

##### PLOT
fig = plt.figure( facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection='3d')
col=ax.scatter(x_pos, y_pos,z_pos , c=peakamp,s=30,  vmin=min(peakamp), vmax=max(peakamp),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col)

#ax.scatter(decay[0], decay[1], decay[2] , c='b',s=30,  marker='o')#cmap=cm.gnuplot2_r) # decay point if needed


ax.set_xlabel('positions along NS (m)')
ax.set_ylabel('positions along EW (m)')
ax.set_zlabel('positions Up (m)')
ax.set_title('peak amplitude')
ax.view_init(elev=21., azim=90)#-146.)




fig1 = plt.figure( facecolor='w', edgecolor='k')
ax1 = fig1.add_subplot(111)

col1=ax1.scatter(x_pos, y_pos,30, c=peakamp,  vmin=min(peakamp), vmax=max(peakamp),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col)

   #plt.scatter(sim_antenna_position[bestsim][:,0]+core_x,sim_antenna_position[bestsim][:,1]+core_y,10, c=sim_tot_power[bestsim,:]*p_ratio[bestsim],vmax=maxp, vmin=0,cmap=cm.gnuplot2_r)


#fig1 = plt.figure( facecolor='w', edgecolor='k')
#ax1 = fig.add_subplot(111, projection='3d')
#col1=ax.scatter(x_pos, y_pos,z_pos , c=peakamp,s=30,  vmin=min(peakamp), vmax=max(peakamp),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
#plt.colorbar(col1)

#ax1.set_xlabel('positions along NS (m)')
#ax1.set_ylabel('positions along EW (m)')
#ax1.set_zlabel('positions Up (m)')
#ax1.set_title('peak amplitude')
#ax1.view_init(elev=21., azim=-146.)


#fig2 = plt.figure( facecolor='w', edgecolor='k')
#ax2 = fig.add_subplot(111, projection='3d')
#col2=ax.scatter(x_pos, y_pos,z_pos , c=peakamp,s=30,  vmin=min(peakamp), vmax=max(peakamp),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
#plt.colorbar(col2)

#ax2.set_xlabel('positions along NS (m)')
#ax2.set_ylabel('positions along EW (m)')
#ax2.set_zlabel('positions Up (m)')
#ax2.set_title('peak amplitude')
#ax2.view_init(elev=21., azim=-146.)


plt.show()


'''

###########################################
##### VOLTAGE
thr_cons= 150 # muV
thr_agg=50 # muV

peaktopeak= np.zeros(number_ant)
peaktopeak2= np.zeros(number_ant)
for i in range(0, number_ant): # loop over all antennas
    try: 
        txt=np.loadtxt(path+ 'out_'+str(i)+'.txt')
        peaktopeak[i]=np.max(txt.T[1]+txt.T[2]) +np.abs(np.min(txt.T[1]+txt.T[2])) # Get Peak to Peak value
        #print peaktopeak[i]
    except:
        peaktopeak[i]=0.
        
        
    ####second set    
    try: 
        txt=np.loadtxt(path2+ 'out_'+str(i)+'.txt')
        peaktopeak2[i]=np.max(txt.T[1]+txt.T[2]) +np.abs(np.min(txt.T[1]+txt.T[2])) # Get Peak to Peak value
        #print peaktopeak[i]
    except:
        peaktopeak2[i]=0.
        
#diff_peaktopeak=peaktopeak-peaktopeak2
    
ind_cons=np.where(peaktopeak>=thr_cons)
ind_agg=np.where(peaktopeak>=thr_agg)
#print ind_cons

print "aggressive: ", len(peaktopeak[ind_agg]), " trigger: ", int(len(peaktopeak[ind_agg])/128)
print "conservative : ", len(peaktopeak[ind_cons]), " trigger: ", int(len(peaktopeak[ind_cons])/128)#, ind_cons


#fig2 = plt.figure()



ax2 = fig.add_subplot(132, projection='3d')
col2=ax2.scatter(x_pos[ind_cons], y_pos[ind_cons],z_pos[ind_cons] , c=peaktopeak[ind_cons],s=30,  vmin=min(peaktopeak), vmax=max(peaktopeak),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col2)

ax2.set_xlabel('positions along NS (m)')
ax2.set_ylabel('positions along EW (m)')
ax2.set_zlabel('positions Up (m)')
ax2.set_title('voltage, conservative, trigger: ' +str(int(len(peaktopeak[ind_cons])/128))+ ', '+ str(len(peaktopeak[ind_cons]))+' ant.')
ax2.view_init(elev=21., azim=-146.)


#fig3 = plt.figure()
ax3 = fig.add_subplot(133, projection='3d')
col3=ax3.scatter(x_pos[ind_agg], y_pos[ind_agg],z_pos[ind_agg] , c=peaktopeak[ind_agg],s=30,  vmin=min(peaktopeak), vmax=max(peaktopeak),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col3)

ax3.set_xlabel('positions along NS (m)')
ax3.set_ylabel('positions along EW (m)')
ax3.set_zlabel('positions Up (m)')
ax3.set_title('voltage, aggressive, trigger: ' +str(int(len(peaktopeak[ind_agg])/128))+ ', '+ str(len(peaktopeak[ind_agg]))+' ant.')
ax3.view_init(elev=21., azim=-146.)

#plt.show()

################ Difference
diff_amp=(peakamp-peakamp2)/max(peakamp2)     # diff refer to max of actual sim
diff_peaktopeak=(peaktopeak-peaktopeak2)/max(peaktopeak2)
print diff_peaktopeak

##### PLOT
figb = plt.figure(figsize=(16, 7), facecolor='w', edgecolor='k')
axb = figb.add_subplot(131, projection='3d')
colb=axb.scatter(x_pos, y_pos,z_pos , c=diff_amp,s=30,  vmin=min(diff_amp), vmax=max(diff_amp),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(colb)

axb.set_xlabel('positions along NS (m)')
axb.set_ylabel('positions along EW (m)')
axb.set_zlabel('positions Up (m)')
axb.set_title('peak amplitude')
axb.view_init(elev=21., azim=-146.)



ax2b = figb.add_subplot(132, projection='3d')
col2b=ax2b.scatter(x_pos[ind_cons], y_pos[ind_cons],z_pos[ind_cons] , c=diff_peaktopeak[ind_cons],s=30,  vmin=min(diff_peaktopeak), vmax=max(diff_peaktopeak),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col2b)

ax2b.set_xlabel('positions along NS (m)')
ax2b.set_ylabel('positions along EW (m)')
ax2b.set_zlabel('positions Up (m)')
ax2b.set_title('voltage, conservative, trigger: ' +str(int(len(peaktopeak[ind_cons])/128))+ ', '+ str(len(peaktopeak[ind_cons]))+' ant.')
ax2b.view_init(elev=21., azim=-146.)


#fig3 = plt.figure()
ax3b = figb.add_subplot(133, projection='3d')
col3b=ax3b.scatter(x_pos[ind_agg], y_pos[ind_agg],z_pos[ind_agg] , c=diff_peaktopeak[ind_agg],s=30,  vmin=min(diff_peaktopeak), vmax=max(diff_peaktopeak),  marker='o', cmap=cm.gnuplot2_r)#cmap=cm.gnuplot2_r)
plt.colorbar(col3b)

ax3b.set_xlabel('positions along NS (m)')
ax3b.set_ylabel('positions along EW (m)')
ax3b.set_zlabel('positions Up (m)')
ax3b.set_title('voltage, aggressive, trigger: ' +str(int(len(peaktopeak[ind_agg])/128))+ ', '+ str(len(peaktopeak[ind_agg]))+' ant.')
ax3b.view_init(elev=21., azim=-146.)

plt.show()
'''









