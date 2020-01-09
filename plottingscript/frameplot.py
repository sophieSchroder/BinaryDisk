#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib import rc 
rc('text',usetex=True)
rc('font',family='serif')

#reading the hdf5 data dumps
frame0 = athdf("BD_comoving.out1.00019.athdf")
frame1 = athdf("BD_nonfollow.out1.00019.athdf")

#the range in R direction we want to plot
rmin = 0
rmax = 128
#cut theta in the middle
thetacut = np.shape(frame0['x2f'])[0]/2

#creat np.meshgrid from the coordinate
#['x1f'] is the face-centered R direction coordinates, ['x3f'] is the face-centered Phi direction coordinates. 
r,phi = np.meshgrid(frame0['x1f'][rmin:rmax], frame1['x3f'])
X = r*np.sin(phi)
Y = r*np.cos(phi)

#Plotting
fig,ax = plt.subplots(1,2)
ax0 = ax[0]
ax1 = ax[1]

ax0.pcolormesh(Y,X,frame0['rho'][:,thetacut,:],vmin=0.01, vmax=1.4) 
ax0.set_aspect(1.0)
ax0.set_title("corotating frame, time="+str(frame0['Time']))
ax0.set_xlabel(r"$r\cos\phi$")
ax0.set_ylabel(r"$r\sin\phi$")

ax1.pcolormesh(Y,X,frame1['rho'][:,thetacut,:],vmin=0.01, vmax=1.4) 
ax1.set_aspect(1.0)
ax1.set_title("non-corotating frame"+str(frame1['Time']))
ax1.set_xlabel(r"$r\cos\phi$")
ax1.set_ylabel(r"$r\sin\phi$")

plt.show()



