#!/usr/bin/env python
import sys
import numpy as np


#INPUT
caseName = "CP1"
nx = 10
ny = 10
nz = 10
Lx = 1.0
Ly = 1.0
Lz = 1.0
nc = nx*ny*nz
shift = Lx/nx
np.random.seed(seed=233423)

# Add some random shifts. 
shiftz = shift*np.random.random_sample((nx*ny,4,nz+1)) - shift/2
shiftx = 0.2*np.random.random_sample(nx-1) -0.1
shifty = 0.2*np.random.random_sample(ny-1) -0.1
#print shiftx

# Write the dimentions. 
dim = open("%s.DIM" % caseName, 'w+')
dim.write("DIMENS")
dim.write('\n')
dim.write('%d %d %d /\n' % (nx, ny, nz))

# Write the COORDS
#print shiftz
#shiftcoord = np.random(nx*ny*nz*8,1)

xcoord = np.linspace(0,Lx,nx+1)
ycoord = np.linspace(0,Ly,ny+1)
zcoord = np.linspace(0,Lz,2)
#print xcoord
#for i in range(0,nx-1):
#    xcoord[i+1] = xcoord[i+1] + shiftx[i,1]
#for i in range(0,ny-1):
#    ycoord[i+1] = ycoord[i+1] + shifty[i,1]
#print xcoord
#print ycoord

grid = open("%s.GRID" % caseName, 'w+')
grid.write("COORD")
grid.write('\n')
for y in np.nditer(ycoord):
    for x in np.nditer(xcoord):
        for z in np.nditer(zcoord):
	    grid.write('%.4f ' % x)
	    grid.write('%.4f ' % y)
	    grid.write('%.4f ' % z)
        grid.write('\n')
grid.write('/\n')

#Write the ZCORN
pillar1 = np.linspace(0,Lz,nz+1)
grid.write("ZCORN")
grid.write('\n')
for k in range(0,nz+1):

    if k == 0:
	grid.write('%d*%.1f \n' % (nx*ny*4, 0.0))
        continue

    if k == nz:
	grid.write('%d*%.1f \n' % (nx*ny*4, Lz))
        continue

    for i in range(0,nx*ny):
        for j in range(0,4):
            top = (pillar1[:k] + shiftz[i][j][:k]).max()
	    bot = (pillar1[k+1:] + shiftz[i][j][k+1:]).min()
            #print top
            #print bot
            val = pillar1[k] + shiftz[i][j][k]
            #print val
            val = min(max(val,top),bot)
            val = min(max(val, 0), Lz)
            #print val
            #val = pillar1[k] + shiftz[i*j][k]
            grid.write('%.4f ' % val)	
	grid.write('\n')
    grid.write('\n')

    for i in range(0,nx*ny):
        for j in range(0,4):
            top = (pillar1[:k] + shiftz[i][j][:k]).max()
	    bot = (pillar1[k+1:] + shiftz[i][j][k+1:]).min()
            #print top
            #print bot
            #print pillar1[k]
            val = pillar1[k] + shiftz[i][j][k]
            val = min(max(val,top),bot)
            val = min(max(val, 0), Lz)
            grid.write('%.4f ' % val)
        grid.write('\n')
    grid.write('\n')
grid.write('/\n')

# Write poro, perm fields to make in runnable in Flow
grid.write("PORO")
grid.write('\n')
grid.write('%d*%.1f / \n' % (nc, 0.3))

grid.write("PERMX")
grid.write('\n')
grid.write('%d*%.1f / \n' % (nc, 100))

grid.write("PERMY")
grid.write('\n')
grid.write('%d*%.1f / \n' % (nc, 100))

grid.write("PERMZ")
grid.write('\n')
grid.write('%d*%.1f / \n ' % (nc, 100))
