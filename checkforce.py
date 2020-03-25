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
#import pandas as pd 
rc('text',usetex=True)
rc('font',family='serif',size=12)


#reading the hdf5 data dumps
frame = athdf("thindisk.out1.00079.athdf")

#the range in R direction we want to plot
rmin = 0
rmax = 384

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x2f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)

rc,phic = np.meshgrid(frame['x1v'][rmin:rmax], frame['x2v'])

#Plotting
fig = plt.figure(figsize=(15,7.5))
spec = gridspec.GridSpec(ncols=2,nrows=1)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,-1])
#cbar = fig.add_subplot(spec[0,-1])

pgas = frame['press'][0,:,rmin:rmax]
gradpgas = np.gradient(pgas, axis=1)
dr = np.gradient(frame['x1v'])
gradpgas /= dr
grav = 0.75/(rc*rc)

im = ax0.pcolormesh(Y,X,gradpgas, cmap='bwr', norm=colors.SymLogNorm(linthresh=1.0e-3))
ax0.set_aspect(1.0)
ax0.set_title("time="+str(frame['Time']))
ax0.set_xlabel(r"$r\cos\phi$")
ax0.set_ylabel(r"$r\sin\phi$")

im2 = ax1.pcolormesh(Y,X,gradpgas/grav, cmap='bwr', norm=colors.SymLogNorm(linthresh=1.0e-3))
ax1.set_aspect(1.0)
ax1.set_title("time="+str(frame['Time']))
ax1.set_xlabel(r"$r\cos\phi$")
ax1.set_ylabel(r"$r\sin\phi$")

#fig.colorbar(im,cax=cbar)
#cbar.set_xlabel(r"$\rho$")
    
plt.savefig('force.png')

