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


#creat movie writter
FFwriter = animation.FFMpegWriter(fps=12, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])

#reading the hdf5 data dumps
frame = athdf("BDstream.out1.00155.athdf")

#the range in R direction we want to plot
rmin = 0
rmax = 256

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x2f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
#X = r*np.sin(phi)
#Y = r*np.cos(phi)

X = np.log10(r)
Y = phi-np.pi

#Plotting
fig = plt.figure(figsize=(7.5,7.5))
spec = gridspec.GridSpec(ncols=30,nrows=1)
ax0 = fig.add_subplot(spec[0,:-1])
cbar = fig.add_subplot(spec[0,-1])

im = ax0.pcolormesh(X,Y,frame['rho'][0,:,rmin:rmax],norm=colors.LogNorm(vmin=0.01, vmax=0.4), cmap='gist_heat') 
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
    


plt.show()