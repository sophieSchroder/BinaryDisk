#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import animation
from matplotlib import rc 
import pandas as pd

rc('text',usetex=True)
rc('font',family='serif',size=12)

GM1 = 0.23076923 
GM2 = 0.76923077
PI = 6.283185307179586/2
Omega = 1.0

#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 192*2
thmin = 0
thmax = 352*2



def plotvalues(csv, GM1, color, ls, label, ax0, ax1, ax2, ax3):

    csv = pd.read_csv(csv)

    rad_avg = csv.rad_avg
    radius = csv.radius
    Dens = csv.surfacedens
    vkep = np.sqrt(GM1/rad_avg)
    AMkep = rad_avg*vkep*Dens
    AM = csv.AM
    MDot = csv.MDot
    AlphaEff = csv.alphaeff



    ax0.plot(radius, Dens,  lw=1.5, c=color, ls=ls, label=label)
    ax0.set_ylabel(r"$\Sigma$")
    ax0.set_xlabel(r"R/L1")
    ax0.set_xlim(np.min(radius), np.max(radius))
    ax0.set_xscale("log")
    ax0.set_ylim(0.001,100.0)
    ax0.set_yscale("log")
    ax0.legend()

    ax1.plot(radius, AM, lw=1.5, c=color, ls=ls)
    ax1.plot(radius, AMkep, lw=1.0, c=color, ls=':')
    ax1.set_ylabel(r"AM")
    ax1.set_xlabel(r"R/L1")
    ax1.set_xlim(np.min(radius), np.max(radius))
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax2.plot(radius, MDot, lw=1.5, c=color, ls=ls)
    ax2.set_ylabel(r"$\dot{M}$")
    ax2.set_xlabel(r"R/L1")
    ax2.set_xlim(np.min(radius), np.max(radius))
    ax2.set_xscale("log")
    ax2.set_ylim(0,0.0025)

    ax3.plot(radius, AlphaEff, lw=1.5, c=color, ls=ls)
    ax3.set_ylabel(r"$\alpha_{\rm eff}$")
    ax3.set_xlabel(r"R/L1")
    ax3.set_xlim(np.min(radius), np.max(radius))
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylim(0.008,1.5)


#Plotting
fig = plt.figure(figsize=(12,8.5))
fig.suptitle(r"t=125-150~for~non-moving~stream, t=97-117~for~moving~stream", fontsize=15)
spec = gridspec.GridSpec(ncols=2,nrows=2)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])
ax2 = fig.add_subplot(spec[1,0])
ax3 = fig.add_subplot(spec[1,1])

#plotvalues('./q5circ/q5circ.csv', 0.16666667, 'tab:purple', '-', r'ecc=0.0, q=5', ax0, ax1, ax2, ax3)
#plotvalues('./q5ecc02/q5e02.csv', 0.16666667, 'tab:cyan', '-', r'ecc=0.2, q=5', ax0, ax1, ax2, ax3)
plotvalues('./q5test.csv', 0.16666667, 'tab:green', '-', r'ecc=0.2, q=5,~moving~stream', ax0, ax1, ax2, ax3)
#plotvalues('./movingstream_ecc035/q5e035_movstr.csv', 0.16666667, 'tab:blue', '-', r'ecc=0.35, q=5,~moving~stream', ax0, ax1, ax2, ax3)


#plt.savefig("compare_boundary_diskprofile.png")
plt.show()

