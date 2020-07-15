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
rc('font',family='serif',size=15)

GM1 = 0.7692307692307692 
PI = 6.283185307179586/2

#creat movie writter

#reading the hdf5 data dumps
frame = athdf("BDstream.out1.00080.athdf")

#the range in R direction we want to plot
rmin = 0
rmax = 384*2

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x2f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)

X_ = np.log10(r)
Y_ = phi

#Plotting
fig = plt.figure(figsize=(16,7.5))
spec = gridspec.GridSpec(ncols=60,nrows=1)
ax0 = fig.add_subplot(spec[0,0:27])
ax1 = fig.add_subplot(spec[0,30:57])
cbar = fig.add_subplot(spec[0,-1])

rho = frame['rho'][0,:,rmin:rmax]

im = ax0.pcolormesh(Y,X,rho, cmap='gist_heat',norm=colors.LogNorm(vmin=0.01, vmax=3.0)) 
ax0.set_xlabel(r"$r\sin\phi$")
ax0.set_ylabel(r"$r\cos\phi$")
ax0.set_aspect('auto')
ax0.set_xlim(np.min(X), np.max(X))
ax0.set_ylim(np.min(Y), np.max(Y))


ax1.pcolormesh(X_,Y_,rho, cmap='gist_heat',norm=colors.LogNorm(vmin=0.01, vmax=3.0)) 
ax1.set_xlabel(r"$\rm Log~R$")
ax1.set_ylabel(r"$\phi$")
ax1.set_xlim(np.min(X_), np.max(X_))
ax1.set_ylim(np.min(Y_), np.max(Y_))

fig.colorbar(im,cax=cbar)
fig.suptitle(r"$\gamma=1.1$"+",~time="+str(frame['Time']), fontsize=20)
cbar.set_ylabel(r"$\rm Log\Sigma$")
    
plt.savefig('density.png')

#plt.show()
