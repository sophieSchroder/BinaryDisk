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

#Plotting
fig = plt.figure(figsize=(20,7.5))
spec = gridspec.GridSpec(ncols=2,nrows=1)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])



def plotam_nonrotate(csvname, ax0):

    data = pd.read_csv(csvname)

    rad_avg = data.rad_avg
    dAM = data.dAM
    AMMdot = data.AMMdot
    AMTH_r = data.AMTH_r
    AMTH_nr = data.AMTH_nr
    TR_r = data.TR_r
    TR_nr = data.TR_nr
    TRdiff = data.TRdiff

    ax0.plot(rad_avg, dAM, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH_nr, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR_nr, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH_nr+TR_nr), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH_nr+TR_nr), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")

    ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    ax0.legend(loc='best', frameon=False)
    ax0.set_title("non-rotating frame")


def plotam_rotate(csvname, ax0):
    data = pd.read_csv(csvname)

    rad_avg = data.rad_avg
    dAM = data.dAM
    AMMdot = data.AMMdot
    AMTH_r = data.AMTH_r
    AMTH_nr = data.AMTH_nr
    TR_r = data.TR_r
    TR_nr = data.TR_nr
    TRdiff = data.TRdiff

    ax0.plot(rad_avg, dAM, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH_r, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR_r, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH_r+TR_r), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH_r+TR_r), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")


    ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    #ax0.set_ylim(-1.e-6, 1.5e-5)
    ax0.legend(loc='best', frameon=False)
    ax0.set_title("co-rotating frame")


plotam_nonrotate("diskscale01_amdt.csv", ax0)
plotam_rotate("diskscale01_amdt.csv", ax1)

plt.show()
