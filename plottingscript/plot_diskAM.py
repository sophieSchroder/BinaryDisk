#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf #need athena_read.py, it's in ~/code/vis/python
import numpy as np
import matplotlib.pyplot as plt
from athena_read import athdf
import matplotlib.gridspec as gridspec
from matplotlib import rc 
import pandas as pd

rc('text',usetex=True)
rc('font',family='serif',size=12)

##################
#XS: there are three frames: 
#frame1. the frame in which both objects are rotating around the center of mass. 
#frame2. the frame in which the primary is fixed, the secondary moves around primary, f_ext only has gravity terms
#frame3. the frame in which both objects are fixed, f_ext also includes inertial forces
#We call frame2 "non-rotating" frame, frame3 "co-rotating frame"
#This script plots AM budget in both frames.
#################

#read in initial conditions
GM1 = 0.16666667
GM2 = 0.83333333
PI = 6.283185307179586/2
Omega = 1.0

#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 192*2
thmin = 0
thmax = 352*2

# #Plotting
# fig = plt.figure(figsize=(20,7.5))
# spec = gridspec.GridSpec(ncols=2,nrows=1)
# ax0 = fig.add_subplot(spec[0,0])
# ax1 = fig.add_subplot(spec[0,1])

fig, axs = plt.subplots(4,1, sharex=True, figsize=(8,10), \
                        gridspec_kw={'height_ratios': [2.2, 1, 1, 1]})
fig.subplots_adjust(hspace=0)
ax0 = axs[0]
ax1  =axs[1]
ax2 = axs[2]
ax3 = axs[3]

def plotam_nonrotate(csvname, ax0):

    data = pd.read_csv(csvname)

    rad_avg = data.rad_avg/np.max(data.rad_avg)
    dAM_r = data.dAM_r
    dAM_nr = data.dAM_nr
    AMMdot = data.AMMdot
    AMTH_r = data.AMTH_r
    AMTH_nr = data.AMTH_nr
    TR_r = data.TR_r
    TR_nr = data.TR_nr
    TRdiff = data.TRdiff

    ax0.plot(rad_avg, dAM_nr, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH_nr, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR_nr, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH_nr+TR_nr), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH_nr+TR_nr), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")

    #ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    ax0.legend(loc='lower left', frameon=False)
    ax0.text(0.05, 0.9, "non-rotating frame", transform=ax0.transAxes)
    ax0.set_xscale('log')
    ax0.set_ylabel(r"$\dot{AM}$")
    ax0.set_ylim(-7.0e-6, 7.0e-6)


def plotam_rotate(csvname, ax0):
    data = pd.read_csv(csvname)

    rad_avg = data.rad_avg/np.max(data.rad_avg)
    dAM_r = data.dAM_r
    dAM_nr = data.dAM_nr
    AMMdot = data.AMMdot
    AMTH_r = data.AMTH_r
    AMTH_nr = data.AMTH_nr
    TR_r = data.TR_r
    TR_nr = data.TR_nr
    TRdiff = data.TRdiff

    ax0.plot(rad_avg, dAM_r, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH_r, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR_r, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH_r+TR_r), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH_r+TR_r), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")


    #ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    #x0.set_ylim(-1.e-6, 1.5e-5)
    ax0.legend(loc='lower left', frameon=False)
    ax0.text(0.05, 0.9, "co-rotating frame", transform=ax0.transAxes)
    ax0.set_xscale('log')
    ax0.set_ylabel(r"$\dot{AM}$")
    ax0.set_ylim(-7.0e-6, 7.0e-6)


def plotvalues(csv, GM1, color, ls, label, ax0, ax2, ax3):

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
    #ax0.legend()

    ax2.plot(radius, MDot, lw=1.5, c=color, ls=ls)
    ax2.set_ylabel(r"$\dot{M}$")
    ax2.set_xlabel(r"R/L1")
    ax2.set_xlim(np.min(radius), np.max(radius))
    ax2.set_xscale("log")
    #ax2.set_ylim(0,0.0025)

    ax3.plot(radius, AlphaEff, lw=1.5, c=color, ls=ls)
    ax3.set_ylabel(r"$\alpha_{\rm eff}$")
    ax3.set_xlabel(r"R/L1")
    ax3.set_xlim(np.min(radius), np.max(radius))
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylim(0.008,1.5)

#plotam_nonrotate("diskscale01_amdt.csv", ax0)
plotam_rotate("diskscale01_amdt.csv", ax0)
plotvalues('q5localiso_75100.csv', 0.16666667, \
           'k', '-', r'', ax1, ax2, ax3)

plt.savefig("AMdisk_h01.png")
plt.show()
