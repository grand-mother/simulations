
# this python script should mimic the functionaliyt of processZhairesShowers.m
# this means it reads in the original ZHAires_output timefresnel_root.dat + inp-file and splits it into single antenna files a_#.dat. It creates also antpos.dat with all antenna positions

import sys
from sys import argv

import numpy as np
import pylab as pl

#import os
#import re


def split(wkdir, steerfile_sim):
    pl.ion()
    DISPLAY = 0

    fname = wkdir+'timefresnel-root.dat' #ZHAires_output

    print 'file ' +fname+ ' gets split up into single antenna file'

    # First load file
    a = np.loadtxt(fname, dtype='float', comments='#')
    if len(a[0]) > 5: # old output format containing vector field etc
        shId = a[:,0]
        antId =  a[:,1]
        x =   a[:,2]  #X = S->N
        y =   a[:,3]  #Y = E->W
        z =   a[:,4]  # Z = Up
        t =   a[:,5]  #ns 
        Ex =  a[:,11]  # V/m 
        Ey =  a[:,12]
        print 'Ez and zpos flipped to correct for coordinate system used for simulations'
        Ez = -1.* a[:,13]

    if len(a[0]) == 5: # new output format
        antId =  a[:,0]
        t =   a[:,1]  #ns 
        Ex =  a[:,2]  # V/m 
        Ey =  a[:,3]
        print 'Ez and zpos flipped to correct for coordinate system used for simulations'
        Ez = -1.* a[:,4] # if RASPASS upgoing needed


        

########## file for antenna positions -- flipping z axis not needed anzmore since read in from inp file
    file_antpos= wkdir+"antpos.dat"
    FILE2 = open(file_antpos, "w" )
    
    datafile = file(steerfile_sim)
    npos=0
    for line in datafile:
        #if 'RASPASSHeight' in line:
        if 'AddAntenna' in line:
          npos+=1
          x_pos = float(line.split(' ',-1)[1])
          y_pos = float(line.split(' ',-1)[2])
          z_pos = float(line.split(' ',-1)[3]) 
          #print x_pos, y_pos, z_pos
          try: 
            print >>FILE2,"%.2f	%.2f	%.2f" % (x_pos , y_pos,z_pos )

          except IndexError: # catch the error
            continue 
    FILE2.close()
    
    
    print "time in ns, Efield in muV/m" 
    for i in range(0, npos):
      sel = np.where(antId == i+1)[0]

######## electric field trace
      ti = t[sel] #*1e-3  # [ns] 
      Exi = Ex[sel]*1e6   # [muV/m]
      Eyi = Ey[sel]*1e6   # [muV/m]
      Ezi = Ez[sel]*1e6   # [muV/m]

      filename = wkdir+"a"+str(i)+".trace"
      alld = np.transpose([ti,Exi,Eyi,Ezi])
      np.savetxt(filename,alld,fmt='%.1f  %.2e  %.2e  %.2e')  
    
    
      #if DISPLAY:
        #ti0 = ti-ti[0]
        #pl.figure(1)
        #pl.plot(ti0,Exi,label='Ex (North-South)')
        #pl.plot(ti0,Eyi,label='Ey (East-West)')
        #pl.plot(ti0,Ezi,label='Ez (Up-Down)')
        #pl.xlabel('Time [mus]')
        #pl.ylabel('Amplitude [muV/m]') 
        #pl.grid(True)
        #pl.legend(loc='lower right');
        #pl.show()
        #raw_input()

    print 'The end'

if __name__ == '__main__':


    wkdir = sys.argv[1] # path where the simulation file is
    steerfile_sim= sys.argv[2] # path to inp-file

    split(wkdir,steerfile_sim)