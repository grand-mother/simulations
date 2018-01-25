#compute the zenith and azimuth angles of the vector normal to the ground at each antenna position
#(antenna xy coordinates are store in a presumed txt file)
#it uses topography (it needs to be linked to the topography framework of valentin)
#all distances in metres

#distance to calculate the local slope for each antenna
deltad=100 #m

#input file with antenna x and y coordinates
antposfile='.txt'
text=np.loadtxt(antposfile)
ant=text[:,0]
xant=text[:,1]
yant=text[:,2]
zant=np.zeros(len(xant))

#loop on all antennas to calculate z, alpha (zenith) and beta (azimuth) angles
for i in range(len(zant)):
    zant[i]=topo.ground_altitude(xant[i],yant[i]):q
    zdeltadx=topo.ground_altitude(xant[i]+deltad,yant[i])
    zdeltady=topo.ground_altitude(xant[i],yant[i]+deltad)
    u=[deltad, 0, zdeltadx-zant[i]]
    v=[0, deltad, zdeltady-zant[i]]
    w=np.cross(u,v)
    w=w/np.linalg.norm(w)
    alpha=np.arccos(w[2]) #zenithal in deg
    beta=np.arccos(w[0]/np.sin(alpha))*180/np.pi #azimuthal in deg
    alpha=alpha*180/np.pi

#output file with antenna xyz coordinates, alpha (zenith) and beta (azimuth) angles
antposfileout='.txt'
with open(antposfileout,'w') as f:
    for i in range(len(zant)):
        f.write(str(ant)+' '+str(xant)+' '+str(yant)+' '+str(zant)+' '+str(alpha)+' '+str(beta)+'\n')
