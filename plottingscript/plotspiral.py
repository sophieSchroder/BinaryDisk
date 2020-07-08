#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf
import numpy as np
from time import time
import multiprocessing as mp
import matplotlib.pyplot as plt
from athena_read import athdf
import matplotlib.gridspec as gridspec
from matplotlib import animation
from matplotlib import rc 

rc('text',usetex=True)
rc('font',family='serif',size=12)

GM1 = 0.7692307692307692 

#creat movie writter

#reading the hdf5 data dumps
frame = athdf("BDstream.out2.00150.athdf")

#the range in R direction we want to plot
rmin = 0
rmax = 384

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x2f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)

#Plotting
fig = plt.figure(figsize=(7.5,7.5))
spec = gridspec.GridSpec(ncols=30,nrows=1)
ax0 = fig.add_subplot(spec[0,:-1])
cbar = fig.add_subplot(spec[0,-1])

rho = frame['rho'][0,:,rmin:rmax]

im = ax0.pcolormesh(Y,X,rho) 
ax0.set_aspect(1.0)
ax0.set_title("time="+str(frame['Time']))
#ax0.set_xlabel(r"$r\sin\phi$")
#ax0.set_ylabel(r"$r\cos\phi$")
ax0.set_xlabel(r"$\rm Log~R$")
ax0.set_ylabel(r"$\phi$")
ax0.set_aspect('auto')
ax0.set_xlim(np.min(X), np.max(X))
ax0.set_ylim(np.min(Y), np.max(Y))

fig.colorbar(im,cax=cbar)
cbar.set_ylabel(r"$\rm Log\Sigma$")
    
#plt.savefig('spiral.pdf')

plt.show()
