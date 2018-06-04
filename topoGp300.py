#import matplotlib
#matplotlib.use('Agg')
import random
import numpy
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import ndimage as ndimage

sys.path.append('/home/renault/retro/lib/')
sys.path.append('/home/renault/retro/lib/python/')
from grand_tour import Topography
sys.path.append('/home/renault/retro/deps/turtle/src')
import turtle

A = numpy.array([0,0])
B = numpy.array([12000,14000])
u = B-A
traj = A+u*numpy.linalg.norm(u)

#Load the topography
latitude, longitude = 43.5,94.0
print "Loading topography map around (lat,long)=",latitude, longitude,"..."
topo = Topography(latitude=latitude, longitude=longitude,path="share/topography", stack_size=121)
print "Done."


# Plot the topography
step = 1000
step_500 = 500
step_250 = 250
xmin = -8000 #-50000
xmax = +14000 #+50000
ymin = -11000 #-50000
ymax = +15000 #+50000
xt = numpy.linspace(xmin, xmax, 251)
yt = numpy.linspace(ymin, ymax, 251)
zt = numpy.zeros((len(yt), len(xt)))
alphat = numpy.zeros((len(yt), len(xt)))
betat = numpy.zeros((len(yt), len(xt)))
lat1 = 43.612421
long1 = 93.820829
alt1 = 2783
ant0 = topo.lla_to_local(lat1,long1,alt1)
print "First antenna location:",ant0

# Create a mask
xt_grid,yt_grid = numpy.meshgrid(xt,yt)
coef_mask = 1.
Kouest = 6000. #7500.
ylim_est = coef_mask*xt_grid + Kouest
Kest = -8000. #-10000. #-16250
ylim_ouest = coef_mask*xt_grid+Kest
mask = numpy.zeros((len(yt), len(xt)))
mask[(yt_grid<ylim_est) & (yt_grid>ylim_ouest)] = 1

# Read and store the slopes
#zt = map(topo.ground_altitude,xt,yt)
for i, yi in enumerate(yt):
    for j, xj in enumerate(xt):
        zt[i, j] = topo.ground_altitude(xj, yi)
        tmp,alphat[i,j],betat[i,j] = topo.ground_normal(xj,yi,0,1)
#        print j,i,topo.ground_normal(xj,yi,0,1),'deg'
	#lla = topo.local_to_lla([xj, yi,zt[i, j]])  # Get topo
        #latitude, longitude, altitude = lla

alphat2 = numpy.copy(alphat)
ind = numpy.where(alphat2>25)
alphat2[ind]=90
alphat2 = ndimage.filters.gaussian_filter(alphat2,1)
mask[alphat>20]=0

#Build 1km-step array
Nx = 40 #12
Ny = 40 #25
ymin_ant = ymin-step
xmin_ant = xmin-step #0. #xmin
xant = numpy.array([xmin_ant+k*step for k in range(Nx)])
yant = numpy.array([ymin_ant+k*step for k in range(Ny)])
xant_grid,yant_grid = numpy.meshgrid(xant,yant)
xant_i = numpy.reshape(xant_grid,-1)
yant_i = numpy.reshape(yant_grid,-1)

#Build 500m and 250m-step arrays
#ind_mid = numpy.argsort(numpy.sqrt((yant-(ant0[1]-14500))**2+(xant-(ant0[0]-10500))**2))
ind_mid = numpy.argsort(numpy.sqrt((yant_i)**2+(xant_i)**2))
ymid = numpy.mean(yant_i[ind_mid[0]]) + step_500
xmid = numpy.mean(xant_i[ind_mid[0]]) + step_500
xant_500 = numpy.array([xmid+k*step_500 for k in range(-5,6)])
yant_500 = numpy.array([ymid+k*step_500 for k in range(-5,6)])
xant_250 = numpy.array([xmid+k*step_250 for k in range(-3,4)])
yant_250 = numpy.array([ymid+k*step_250 for k in range(-3,4)])
xant_grid,yant_grid = numpy.meshgrid(xant_500,yant_500)
xant_500_i = numpy.reshape(xant_grid,-1)
yant_500_i = numpy.reshape(yant_grid,-1)
xant_grid,yant_grid = numpy.meshgrid(xant_250,yant_250)
xant_250_i = numpy.reshape(xant_grid,-1)
yant_250_i = numpy.reshape(yant_grid,-1)

print 'x,y mid = ',xmid,ymid
print 'initialization finished'

#Initiate criteria
avslope_max = 10.5
ymean_max = 0.
alpha_max = 20.
zmax = 2700.
avslope_ap = 90.
Nant_ap = 0
ymean_ap = ymin
Nant_1km = 175 #162
flag_param = True
while flag_param: #numpy.logical_and(avslope_ap>9.5,ymean_ap<0.):
    #Draw random parameters of the arrays (dx,dy,theta)
    dx = random.uniform(-1.,1.)*step/2.
    dy = random.uniform(-1.,1.)*step/2.
    theta = np.radians(random.uniform(0.,1.)*90)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c,-s), (s, c)))

    #Translate the arrays
    xant = numpy.copy(xant_i)
    yant = numpy.copy(yant_i)
    xant = xant + dx
    yant = yant + dy
    xant_250 = numpy.copy(xant_250_i)
    yant_250 = numpy.copy(yant_250_i)
    xant_500 = numpy.copy(xant_500_i)
    yant_500 = numpy.copy(yant_500_i)
    xant_250 = xant_250 + dx #(dx if (abs(dx)<abs(dx-250)) else (dx-250))
    yant_250 = yant_250 + dy #(dy if (abs(dy)<abs(dy-250)) else (dy-250))
    xant_500 = xant_500 + dx #(dx if (abs(dx)<abs(dx-500)) else (dx-500))
    yant_500 = yant_500 + dy #(dy if (abs(dy)<abs(dy-500)) else (dy-500))

    #Rotate the arrays
    ant = R.dot(numpy.vstack((xant,yant)))
    xant = ant[0,:]
    yant = ant[1,:]
    ant_500 = R.dot(numpy.vstack((xant_500,yant_500)))
    xant_500 = ant_500[0,:]
    yant_500 = ant_500[1,:]
    ant_250 = R.dot(numpy.vstack((xant_250,yant_250)))
    xant_250 = ant_250[0,:]
    yant_250 = ant_250[1,:]
    ant_mid = R.dot(numpy.vstack((xmid,ymid)))
    xmidr = ant_mid[0,:]
    ymidr = ant_mid[1,:]

    #Remove out of bounds antennas in 1km-step array
    ind = numpy.where(numpy.logical_or(yant<(coef_mask*xant+Kest),yant>(coef_mask*xant+Kouest)))
    xant = numpy.delete(xant,ind)
    yant = numpy.delete(yant,ind)
    ind = numpy.where(numpy.logical_or(numpy.logical_or(yant>ymax,yant<ymin),numpy.logical_or(xant>xmax,xant<xmin)))
    xant = numpy.delete(xant,ind)
    yant = numpy.delete(yant,ind)
    ind = numpy.where(numpy.logical_or(numpy.logical_or(yant>ymax,yant<ymin),numpy.logical_or(xant>xmax,xant<xmin)))
    xant = numpy.delete(xant,ind)
    yant = numpy.delete(yant,ind)

    # Filter antennas to avoid isolated and/or too inclined ones
    ix=0
    flag = True
    while flag:
        dist = numpy.sqrt((xant-xant[ix])**2+(yant-yant[ix])**2)
        Ndist = numpy.size(numpy.where(dist<=1.5*step))-1
        Ndist_close = numpy.size(numpy.where(numpy.sqrt((xant-xant[ix])**2+(yant-yant[ix])**2)<=step))
        tmp,alpha_ant,tmp=topo.ground_normal(xant[ix],yant[ix],0,1)
        zant = topo.ground_altitude(xant[ix], yant[ix])
        if (alpha_ant>alpha_max and Ndist<5) or numpy.logical_and(Ndist<3,Ndist_close<2) or zant>zmax: # or alpha_ant>alpha_max+10.:  #7 #4
            yant = numpy.delete(yant,ix)
            xant = numpy.delete(xant,ix)
            ix=0
            flag = True
        else:
            ix=ix+1
        if ix>len(xant)-1:
            flag = False

    #Remove antenna too far south
    Nant = len(xant)
    if Nant>Nant_1km:
        dist_ant0 = numpy.sqrt((xant-ant0[0])**2+(yant-ant0[1])**2)
        ind_dist_ant0 = numpy.argsort(dist_ant0)
        xant = numpy.delete(xant,ind_dist_ant0[-(Nant-Nant_1km):])
        yant = numpy.delete(yant,ind_dist_ant0[-(Nant-Nant_1km):])

    #Compute 1km-step array average slope and easting
    avslope = 0
    zant = numpy.zeros((len(xant)))
    ialpha = 0
    for ix in range(len(xant)):
        tmp,alpha_ant,tmp=topo.ground_normal(xant[ix],yant[ix],0,1)
        zant[ix] = topo.ground_altitude(xant[ix], yant[ix])
        avslope = avslope + alpha_ant
        if alpha_ant>alpha_max:
            ialpha = ialpha+1
    ymean = numpy.mean(yant)
    avslope = avslope / Nant

    # Check slopes and altitudes of 250m- and 500m-step arrays
    flag_500 = True
    alpha_500 = numpy.zeros((len(xant_500)))
    zant_500 = numpy.zeros((len(xant_500)))
    for ix in range(len(xant_500)):
        tmp,alpha_500[ix],tmp = topo.ground_normal(xant_500[ix],yant_500[ix],0,1)
        zant_500[ix] = topo.ground_altitude(xant_500[ix], yant_500[ix])
    alpha_test_500 = numpy.amax(alpha_500)
    z_test_500 = numpy.amax(zant_500)
    if alpha_test_500>alpha_max+10. or z_test_500>zmax:
        flag_500 = False
    if flag_500:
        flag_250 = True
        alpha_250 = numpy.zeros((len(xant_250)))
        zant_250 = numpy.zeros((len(xant_250)))
        for ix in range(len(xant_250)):
            tmp,alpha_250[ix],tmp = topo.ground_normal(xant_250[ix],yant_250[ix],0,1)
            zant_250[ix] = topo.ground_altitude(xant_250[ix], yant_250[ix])
        alpha_test_250 = numpy.amax(alpha_250)
        z_test_250 = numpy.amax(zant_250)
        if alpha_test_250>alpha_max+10. or z_test_250>zmax:
            flag_250 = False
    else:
        flag_250 = False

    # Check centering of 250m- and 500m-step arrays
    if flag_500 and flag_250:
        xmax_500 = numpy.amax(xant_500)
        xmin_500 = numpy.amin(xant_500)
        ind_ymid = abs(yant-ymid)<=2500.
        if numpy.size(numpy.where(xant[ind_ymid]>=xmax_500))>1 and numpy.size(numpy.where(xant[ind_ymid]<=xmin_500))>1:
            flag_center = True
        else:
            flag_center = False

    #Check if this array is worth keeping or not
    if avslope<avslope_ap and Nant>=Nant_1km and flag_500 and flag_250 and flag_center and ialpha<50:
        avslope_ap = avslope
        ymean_ap = ymean
        Nant_ap = Nant
        xant_ap = xant
        yant_ap = yant
        print 'Nant_ap =',Nant_ap,'avslope= ',avslope_ap,'ymean_ap= ',ymean_ap
        print 'alpha_max_250 = ',alpha_test_250,' alpha_max_500 = ',alpha_test_500
        print '1000m-step array :',numpy.shape(xant)[0],' antennas, 500m-step array :',numpy.shape(xant_500)[0],' antennas, 250m-step array:',numpy.shape(xant_250)[0],' antennas.'

        if avslope_ap<avslope_max and ymean_ap>ymean_max and numpy.shape(xant)[0]==175 and numpy.shape(xant_500)[0]==121 and numpy.shape(xant_250)[0]==49:
            print 'all tests successful'
            flag_param = False

print 'Done scanning the array'
print 'Final values ::: Nant =',Nant_ap,' avslope =',avslope_ap,'ymean = ',ymean_ap
print 'dx = ',dx,' dy = ',dy,' theta = ',theta
print ialpha,'one or more antennas with slope>',alpha_max,' deg'

xant = xant_ap
yant = yant_ap

################################
################################
# WARNING : for very horizontal showers, it may run forever because of the shadowing. Introduce a test to prevent it.
################################
################################

#remove common antennas in 250m, 500m-step and 1km-step arrays
ind_del = []
for i250,x250 in enumerate(xant_250):
    for i1000,x1000 in enumerate(xant):
        if x1000==x250 and yant[i1000]==yant_250[i250]:
            ind_del.append(i250)
xant_250 = numpy.delete(xant_250,ind_del)
yant_250 = numpy.delete(yant_250,ind_del)
zant_250 = numpy.delete(zant_250,ind_del)

ind_del = []
for i500,x500 in enumerate(xant_500):
    for i1000,x1000 in enumerate(xant):
        if x1000==x500 and yant[i1000]==yant_500[i500]:
            ind_del.append(i500)
xant_500 = numpy.delete(xant_500,ind_del)
yant_500 = numpy.delete(yant_500,ind_del)
zant_500 = numpy.delete(zant_500,ind_del)

ind_del = []
for i250,x250 in enumerate(xant_250):
    for i500,x500 in enumerate(xant_500):
        if x500==x250 and yant_500[i500]==yant_250[i250]:
            ind_del.append(i250)
xant_250 = numpy.delete(xant_250,ind_del)
yant_250 = numpy.delete(yant_250,ind_del)
zant_250 = numpy.delete(zant_250,ind_del)

#X = SN
#Y = EW
ylim_ouest = xt*coef_mask + Kouest
ylim_est = xt*coef_mask + Kest

print '1000m-step array :',numpy.shape(xant)[0],' antennas, 500m-step array :',numpy.shape(xant_500)[0],' antennas, 250m-step array:',numpy.shape(xant_250)[0],' antennas.'

ant1000 = numpy.vstack((xant,yant,zant))
ant250 = numpy.vstack((xant_250,yant_250,zant_250))
ant500 = numpy.vstack((xant_500,yant_500,zant_500))
antpos = numpy.hstack((ant1000,ant500,ant250)).T
antcoord = numpy.zeros(numpy.shape(antpos))
for iant in range(numpy.shape(antpos)[0]):
    antcoord[iant,:] = topo.local_to_lla(antpos[iant,:])

# Store the results
ant_file = '/home/renault/GP300_antpos.txt'
numpy.savetxt(ant_file,antpos,fmt='%s     %s     %s', header="X [m]   Y [m]   Z [m]")
ant_file = '/home/renault/GP300_antcoord.txt'
numpy.savetxt(ant_file,antcoord,fmt='%s     %s     %s', header="lat [m]   lon [m]   Z [m]")
ant_file = '/home/renault/GP300_array_parameters.txt'
stream = [
"xmin = {:.2f} [m]".format(xmin_ant),
"ymin = {:.2f} [m]".format(ymin_ant),
"xmid = {:.2f} [m]".format(xmid),
"ymid = {:.2f} [m]".format(ymid),
"dx = {:.2f} [m]".format(dx),
"dy = {:.2f} [m]".format(dy),
"theta = {:.2f} [deg]".format(numpy.degrees(theta)),
"average_slope = {:.2f} [deg]".format(avslope_ap),
"average_y = {:.2f} [m]".format(ymean_ap)
]
"\n".join(stream)
param_file = open(ant_file,"w+")
param_file.write(str(stream))
param_file.close()

#Plot the results over the altitude and the slope
plt.figure(1)
ax = plt.gca()
plt.pcolor(xt, yt, zt, cmap="terrain", alpha=0.75)
plt.xlabel("Northing (m)")
plt.ylabel("Westing (m)")
plt.colorbar()
plt.plot(ant0[0],ant0[1],'ks',markersize=4,linewidth=2)
plt.plot(xmidr,ymidr,'ks',markersize=4,linewidth=2)
plt.plot(xant_250,yant_250,'ro',markersize=4,markeredgecolor='r')
plt.plot(xant_500,yant_500,'bo',markersize=4,markeredgecolor='b')
plt.plot(xant,yant,'go',markersize=4,markeredgecolor='g')
plt.plot(xt,ylim_ouest,'k--')
plt.plot(xt,ylim_est,'k--')
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
figname = '/home/renault/array_on_altitude_map.png'
plt.savefig(figname,dpi=750)

plt.figure(2)
ax = plt.gca()
plt.pcolor(xt, yt, alphat2, cmap="terrain", alpha=0.75)
plt.xlabel("Northing (m)")
plt.ylabel("Westing (m)")
plt.colorbar()
plt.plot(ant0[0],ant0[1],'ks',markersize=4,linewidth=2)
plt.plot(xmidr,ymidr,'ks',markersize=4,linewidth=2)
plt.plot(xant_250,yant_250,'ro',markersize=4,markeredgecolor='r')
plt.plot(xant_500,yant_500,'bo',markersize=4,markeredgecolor='b')
plt.plot(xant,yant,'go',markersize=4,markeredgecolor='g')
plt.plot(xt,ylim_ouest,'k--')
plt.plot(xt,ylim_est,'k--')
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
figname = '/home/renault/array_on_slope_map.png'
plt.savefig(figname,dpi=750)

'''
plt.figure(3)
ax = plt.gca()
plt.pcolor(xt, yt, mask, cmap="terrain", alpha=0.75)
plt.xlabel("Northing (m)")
plt.ylabel("Westing (m)")
plt.colorbar()
#plt.plot(ant0[0],ant0[1],'ko')
plt.plot(xant,yant,'ko',markersize=2)
plt.plot(xt,ylim_ouest,'r')
plt.plot(xt,ylim_est,'r')
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
'''

setsA = [[-4000,-6000],[-6000,6000],[3000,20000]] #points de debut
setsB = [[12000,14000],[12000,-9000],[16000,-2000]] #points de fin
col = ["k-","r-","b-"]
for i in range(len(setsA)):
  A = numpy.array(setsA[i])
  B = numpy.array(setsB[i])
  u = B-A

#  plt.figure(1)
#  plt.plot([A[0],B[0]],[A[1],B[1]],col[i])
#  plt.figure(2)
#  plt.plot([A[0],B[0]],[A[1],B[1]],col[i])

  '''
  normu = numpy.linalg.norm(u)
  u = u/normu
  traj = numpy.linspace(0,normu,200)
  xtraj = numpy.array(A[0]+u[0]*traj)
  ytraj = numpy.array(A[1]+u[1]*traj)
  zt = numpy.zeros(len(traj))
  alphat = numpy.zeros(len(traj))
  betat = numpy.zeros(len(traj))
  for j in range(len(traj)):
    zt[j] = topo.ground_altitude(xtraj[j], ytraj[j])
    tmp,alphat[j],betat[j] = topo.ground_normal(xtraj[j],ytraj[j],0,1)

  plt.figure()
  #plt.plot(traj/1e3,alphat,col[i])
  plt.plot(traj/1e3,zt,col[i],linestyle='--')
  plt.grid(True)
  plt.xlabel("Distance along track (km)")
  plt.ylabel("Altitude (m)")
  '''

plt.show()
