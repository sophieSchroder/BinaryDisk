#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf
import numpy as np
from time import time
import multiprocessing as mp
import matplotlib.pyplot as plt
from athena_read import hst
import matplotlib.gridspec as gridspec
from matplotlib import animation
from matplotlib import rc 

rc('text',usetex=True)
rc('font',family='serif',size=12)


#Plotting
fig = plt.figure(figsize=(7.5,10))
spec = gridspec.GridSpec(ncols=1,nrows=3)
ax0 = fig.add_subplot(spec[0,:])
ax1 = fig.add_subplot(spec[1,:])
ax2 = fig.add_subplot(spec[2,:])

def plotmdot(hstname, ax):

    #reading history file
    hstdata = hst(hstname)

    massinflow = hstdata['massfluxix1']-hstdata['massfluxox1']
    cumsum = np.cumsum(massinflow)
    dt = 0.1
    dv = (0.6/128)*(2*np.pi/256)
    time = hstdata['time'][1:]
    mass = hstdata['mass'][1:]

    ax.plot(time, -hstdata['massfluxix1'][1:], label=r'$\dot{M}_{\rm inner}$')
    ax.plot(time, -hstdata['massfluxox1'][1:], label=r'$\dot{M}_{\rm outer}$')
    ax.legend(loc='lower right')
    ax.set_ylabel(r'$\dot{M}$')
    ax.set_xlim(np.min(time), np.max(time))

    axx = ax.twinx()
    axx.set_ylabel(r"$M(t)$", color="tab:red")
    axx.plot(time, hstdata['mass'][1:], color='tab:red', ls='--')
    axx.tick_params(axis='y', labelcolor="tab:red")
    axx.set_ylim(0.0,0.6)


plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/highres_adiabatic/stampede/gamma11_hlle_new/BDstream.hst", ax0)
plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc_fixstream/BDstream.hst", ax1)
plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc035_fixstream/BDstream.hst", ax2)

fig.suptitle(r"$\gamma=1.1$"+r",~fix stream")
ax0.set_title(r"ecc=0.0")
ax1.set_title(r"ecc=0.2")
ax2.set_title(r"ecc=0.35")
ax2.set_xlabel(r"$t/t_{0}$")
#plt.show()
plt.savefig("g11eccpotential_mdot.png")

'''
plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/highres_adiabatic/stampede/gamma11_hlle_new/BDstream.hst", ax0)
plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc_hr/BDstream.hst", ax1)
plotmdot("/LyraShared/xh2pm/BinaryDisk_prep/inflow/sinstream/g11ecc035_sinstream/BDstream.hst", ax2)

fig.suptitle(r"$\gamma=1.1$"+r",~sin stream")
ax0.set_title(r"ecc=0.0")
ax1.set_title(r"ecc=0.2")
ax2.set_title(r"ecc=0.35")
ax2.set_xlabel(r"$t/t_{0}$")
#plt.show()
plt.savefig("g11eccpotential_sinstream_mdot.png")
'''
