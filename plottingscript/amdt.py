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
#import pandas as pd

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
GM1 = 0.7692307692307692 
GM2 = 0.2307692307692307
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


#reading angular momentum
def readam(i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
        name2 = 'BDstream.out2.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']

    AM = np.sum(frame2['user_out_var6'][0,:,:], axis=0)
    AM_net = np.sum(frame2['user_out_var11'][0,:,:], axis=0)
    AM_nc = np.sum(frame2['user_out_var34'][0,:,:], axis=0)
        
    return [Radius, AM, time, AM_net, AM_nc]

def intammdot(i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
        name2 = 'BDstream.out2.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)

    time = frame['Time']
    
    Radius = frame['x1v']
    dR = np.gradient(Radius)
    AMMdot = np.sum(frame2['user_out_var16'][0,:,:], axis=0)

    return [Radius,AMMdot]


def intamfh(i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
        name2 = 'BDstream.out2.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    dR = np.gradient(Radius)
    amfh = np.sum(frame2['user_out_var17'][0,:,:], axis=0)
    amth_net = np.sum(frame2['user_out_var18'][0,:,:], axis=0)
    amth_nc = np.sum(frame2['user_out_var37'][0,:,:], axis=0)

    return [Radius, amfh, amth_net, amth_nc]



def inttorque(i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
        name2 = 'BDstream.out2.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'
        name2 = 'BDstream.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)

    time = frame['Time']
    
    Radius = frame['x1v']
    dR = np.gradient(Radius)
    trr = np.sum(frame2['user_out_var19'][0,:,:], axis=0)
    trr2 = np.sum(frame2['user_out_var20'][0,:,:], axis=0)
        
    return [Radius, trr,  trr2]


def plotam_nonrotate(start, end, ax0):

    rad_avg = readam(end)[0]
    dAM = readam(end)[4]-readam(start)[4]
    AMMdot = intammdot(end)[1] - intammdot(start)[1]
    AMTH_r = intamfh(end)[1] - intamfh(start)[1]
    AMTH_nr = intamfh(end)[3] - intamfh(start)[3]
    TR_nr = inttorque(end)[2] - inttorque(start)[2]
    TR_r  = inttorque(end)[1] - inttorque(start)[1]
    TRdiff = TR_r - TR_nr
    dtime = readam(end)[2]-readam(start)[2]

    dAM /= dtime
    AMMdot /= dtime
    AMTH_r /= dtime
    AMTH_nr /= dtime
    TR_r /= dtime
    TR_nr /= dtime
    TRdiff /= dtime

    ax0.plot(rad_avg, dAM, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH_nr, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR_nr, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH_nr+TR_nr), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH_nr+TR_nr), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")

    ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    ax0.set_ylim(-1.5e-5, 1.5e-5)
    ax0.legend(loc='best', frameon=False)
    ax0.set_title("non-rotating frame")


def plotam_rotate(start, end, ax0):

    rad_avg = readam(end)[0]
    dAM = readam(end)[1]-readam(start)[1]
    AMMdot = intammdot(end)[1] - intammdot(start)[1]
    AMTH = intamfh(end)[1] - intamfh(start)[1]
    TR = inttorque(end)[1] - inttorque(start)[1]
    dtime = readam(end)[2]-readam(start)[2]

    dAM /= dtime
    AMMdot /= dtime
    AMTH /= dtime
    TR /= dtime

    ax0.plot(rad_avg, dAM, label=r"$AM_{t}(R)$", lw=2.0, c='k')
    ax0.plot(rad_avg, AMMdot, lw=2.0, label=r"$AM_{\dot{M}}(R)$", c='y', ls='-.')
    ax0.plot(rad_avg, AMTH, lw=2.0, label=r"$AM_{TH}(R)$", c='deepskyblue')
    ax0.plot(rad_avg, TR, lw=2.0, label=r"$T(R)$", ls='--', c='blue')
    ax0.plot(rad_avg, -(AMTH+TR), lw=2.0, ls=':', label=r'dissipation', c='k')
    ax0.plot(rad_avg, (AMMdot+AMTH+TR), ls='--', lw=2.0, c='r', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")


    ax0.set_xlim(np.min(rad_avg), np.max(rad_avg))
    ax0.set_xlabel(r"R")
    ax0.set_ylim(-1.5e-5, 1.5e-5)
    ax0.legend(loc='best', frameon=False)
    ax0.set_title("co-rotating frame")


plotam_nonrotate(45,20,ax0)
plotam_rotate(45,20,ax1)
plt.show()
#plt.savefig("AMdt_.png")

