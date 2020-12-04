#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf
import numpy as np
import matplotlib.pyplot as plt
from athena_read import athdf
import matplotlib.gridspec as gridspec
from matplotlib import rc 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

rc('text',usetex=True)
rc('font',family='serif',size=15)

GM1 = 0.7692307692307692 
PI = 6.283185307179586/2

#creat movie writter


#the range in R direction we want to plot
rmin = 0
rmax = 384

#Plotting
fig = plt.figure(figsize=(12,5))
spec = gridspec.GridSpec(ncols=2,nrows=1)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])

def plotdensity(name, ax):
    #reading the hdf5 data dumps
    frame = athdf(name)
    rho = frame['rho'][0,:,rmin:rmax]
    r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
    X = r*np.sin(phi)
    Y = r*np.cos(phi)

    im = ax.pcolormesh(Y,X,rho, norm=colors.LogNorm(vmin=0.001, vmax=5.0)) 
    ax.set_xlabel(r"$r\sin\phi$")
    #ax.set_ylabel(r"$r\cos\phi$")
    ax.set_aspect('auto')
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    #ax.set_title(r"t="+str(frame['Time']))
    ax.set_aspect('equal')
    return im

def plotspiral(name, ax):
    #reading the hdf5 data dumps
    frame = athdf(name)
    rho = frame['rho'][0,:,rmin:rmax]
    r,phi = np.meshgrid(frame['x1f'][rmin:rmax], frame['x2f'])
    X = np.log10(r)
    Y = phi

    ax.pcolormesh(X,Y,rho ,norm=colors.LogNorm(vmin=0.001, vmax=3.0)) 
    ax.set_xlabel(r"$\rm Log~R$")
    ax.set_ylabel(r"$\phi$")
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))

im = plotdensity("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc035_fixstream/BDstream.out1.00150.athdf", ax0)
plotdensity("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc035_sinstream/BDstream.out1.00150.athdf", ax1)

ax0.set_title(r"amp = 0.0")
ax1.set_title(r"amp = 0.2")
ax0.set_ylabel(r"$r\cos\phi$")
fig.suptitle(r"$\gamma=1.1$"+r"~ecc=0.35"+r"~t=150")

axins = inset_axes(ax1, width="5%", height="100%", loc="lower left", \
                   bbox_to_anchor=(1.05,0.,1,1),
                   bbox_transform=ax1.transAxes,
                   borderpad=0.2)

fig.colorbar(im, cax=axins)

axins.set_ylabel(r"$\rm log(\rho/\rho_{0})$")
plt.savefig('g11035eccpotential_densitysnapshot.png')


#plt.show()
