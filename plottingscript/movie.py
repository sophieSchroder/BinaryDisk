#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf #need athena_read.py, it's in ~/code/vis/python
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from athena_read import athdf
import matplotlib.gridspec as gridspec
from matplotlib import animation
from matplotlib import rc 
import pandas as pd 
rc('text',usetex=True)
rc('font',family='serif',size=12)


#creat movie writter
FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])

#reading a hdf5 data dumps to get dimension info, change file name as you want : )
frame = athdf("BDstream.out1.00019.athdf")

#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 256

#creat np.meshgrid from the coordinate, for cylindrical
#['x1f'] is the face-centered R direction coordinates, ['x3f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)

#Plotting
fig = plt.figure(figsize=(7.5,7.5))
spec = gridspec.GridSpec(ncols=30,nrows=1)
ax0 = fig.add_subplot(spec[0,:-1])
cbar = fig.add_subplot(spec[0,-1])

def animated(i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'

    frame = athdf(name)
        
    print("plotting "+name+"...")

    im = ax0.pcolormesh(Y,X,frame['rho'][0,:,rmin:rmax],norm=colors.LogNorm(vmin=0.001, vmax=0.5)) 
    ax0.set_aspect(1.0)
    ax0.set_title("time="+str(frame['Time']))
    ax0.set_xlabel(r"$r\cos\phi$")
    ax0.set_ylabel(r"$r\sin\phi$")

    fig.colorbar(im,cax=cbar)
    cbar.set_xlabel(r"$\rho$")
    
    return [im]


anim = animation.FuncAnimation(fig, animated, frames=100, interval=200, blit=True) #change frames= number of frames
anim.save('BDcyl.mov',writer=FFwriter)