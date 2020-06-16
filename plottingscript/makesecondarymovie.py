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
import pandas as pd 

rc('text',usetex=True)
rc('font',family='serif',size=12)

#creat movie writter
FFwriter = animation.FFMpegWriter(fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])

base_dir = "/Users/sophielundschroder/Dropbox/Data_binDisk/ecc_test/ecc/"

#reading in pm_trackfile using pandas
pmtrackfile = pd.read_csv(base_dir  +"pm_trackfile.dat",delim_whitespace=True)

#reading the hdf5 data dumps

frame = athdf(base_dir +"BDrst_ecc.out1.00405.athdf")

#the range in R direction we want to plot, in index
rmin = 0
rmax = 384

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x2f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)


#Plotting
fig = plt.figure(figsize=(12,7.5))
spec = gridspec.GridSpec(ncols=30,nrows=1)
ax0 = fig.add_subplot(spec[0,:-1])
#cbar = fig.add_subplot(spec[0,-1])


def animated(i):
    ax0.clear()
    
    
    #track file time 
    pmtrack_index = i*5+1
    secondary_x = pmtrackfile.x.iloc[1:i+2]
    secondary_y = pmtrackfile.y.iloc[1:i+2]

    #change file names here accordingly
    hdf5_index = 405+i*5

    if hdf5_index < 10:
        name = base_dir +'BDrst_ecc.out1.0000'+str(hdf5_index)+'.athdf'
    elif hdf5_index  >=100:
        name = base_dir +'BDrst_ecc.out1.00'+str(hdf5_index)+'.athdf'
    else:
        name = base_dir +'BDrst_ecc.out1.000'+str(hdf5_index)+'.athdf'

    frame = athdf(name)
        
    print("plotting "+name+"...")

    im = ax0.pcolormesh(Y,X,frame['rho'][0,:,rmin:rmax],norm=colors.LogNorm(vmin=0.001, vmax=0.2),cmap="magma") 
    ax0.scatter(secondary_x, secondary_y, marker='+',c='k',s=2.0)
    ax0.set_aspect('equal')
    ax0.set_title("time="+str(frame['Time']))
    ax0.set_xlabel(r"$r\cos\phi$")
    ax0.set_ylabel(r"$r\sin\phi$")
    ax0.set_ylim(np.min(Y), np.max(Y))
    ax0.set_xlim(np.min(X),1.3)

    #fig.colorbar(im,cax=cbar)
    #cbar.set_xlabel(r"$\rho$")
    
    return [im]


anim = animation.FuncAnimation(fig, animated, frames=240, interval=200, blit=True) #change frames= number of frames
anim.save('BDrst_ecc02.mov',writer=FFwriter)