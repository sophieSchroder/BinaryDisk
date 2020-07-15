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

GM1 = 0.7692307692307692 
GM2 = 0.2307692307692307
PI = 6.283185307179586/2
Omega = 1.0

#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 192*2
thmin = 0
thmax = 352*2

#read in pmtrack
#secondarytrack = pd.read_csv("pm_trackfile.dat",delim_whitespace=True) 

def density_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    dens_avg = np.array([])
    
    for radius in range(rmin,rmax):
        sigma = frame2['user_out_var24'][0,:,radius]
        sigma_avg =  np.sum(sigma)/(2*PI)
        dens_avg = np.append(dens_avg, sigma_avg)
                           
    return [Radius, dens_avg]

def am_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    am_avg = np.array([])
    am_kep = np.array([])
    
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        sigma = rho
        dphi = np.gradient(frame['x2v'])
        am = frame2['user_out_var25'][0,:,radius]
        am = np.sum(am)/(2*PI)
        am_avg = np.append(am_avg, am)
        am_k = np.sum(np.sqrt(GM1*rad)*sigma*dphi)/(2*PI)
        am_kep = np.append(am_kep, am_k)
        
    return [Radius, am_avg, am_kep]

def mdot_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    mdot_avg = np.array([])

    for radius in range(rmin,rmax):
        mdot =  -frame2['user_out_var21'][0,:,radius]
        mdot = np.sum(mdot)
        mdot_avg = np.append(mdot_avg, mdot)

    return [Radius, mdot_avg]


def alphaterm2_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    mdot_avg = np.array([])

    for radius in range(rmin,rmax):
        mdot =  frame2['user_out_var36'][0,:,radius]
        mdot = np.sum(mdot)
        mdot_avg = np.append(mdot_avg, mdot)

    return [Radius, mdot_avg]


def alphaeff_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    alpha_avg = np.array([])

    #rvk = np.sqrt(GM1*Radius) 
    for radius in range(rmin,rmax):

        alpha_eff = frame2['user_out_var22'][0,:,radius]
        alpha_eff  = np.sum(alpha_eff)
        alpha_avg = np.append(alpha_avg, alpha_eff)

    return [Radius, alpha_avg]


def mach_avg(name, name2):

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    mach_avg = np.array([])

    for radius in range(rmin,rmax):
        mach =  frame2['user_out_var33'][0,:,radius]
        mach = np.sum(mach)
        mach_avg = np.append(mach_avg, mach)

    return [Radius, mach_avg]

def plotvalues(gamma, start, end, dir, color, label, ax0, ax1, ax2, ax3):

    #change file names here accordingly
    if start < 10:
        startframename = dir+'BDstream.out1.0000'+str(start)+'.athdf'
        startframename2 = dir+'BDstream.out2.0000'+str(start)+'.athdf'
    elif start >=100:
        startframename = dir+'BDstream.out1.00'+str(start)+'.athdf'
        startframename2 = dir+'BDstream.out2.00'+str(start)+'.athdf'
    else:
        startframename = dir+'BDstream.out1.000'+str(start)+'.athdf'
        startframename2 = dir+'BDstream.out2.000'+str(start)+'.athdf'

    starttime = athdf(startframename)['Time']

    if end < 10:
        endframename = dir+'BDstream.out1.0000'+str(end)+'.athdf'
        endframename2 = dir+'BDstream.out2.0000'+str(end)+'.athdf'
    elif end >=100:
        endframename = dir+'BDstream.out1.00'+str(end)+'.athdf'
        endframename2 = dir+'BDstream.out2.00'+str(end)+'.athdf'
    else:
        endframename = dir+'BDstream.out1.000'+str(end)+'.athdf'
        endframename2 = dir+'BDstream.out2.000'+str(end)+'.athdf'

    endtime = athdf(endframename)['Time']
    dtime = endtime-starttime

    rad_avg = density_avg(startframename, startframename2)[0]
    Dens = density_avg(endframename,endframename2)[1] - density_avg(startframename, startframename2)[1]
    Dens /= dtime
    AM = am_avg(endframename,endframename2)[1] - am_avg(startframename, startframename2)[1]
    AM /= dtime
    MDot = mdot_avg(endframename,endframename2)[1] - mdot_avg(startframename, startframename2)[1]
    MDot /= dtime
    Alpha = alphaeff_avg(endframename,endframename2)[1] - alphaeff_avg(startframename, startframename2)[1]
    Alpha /= dtime
    Alpha_term2 = alphaterm2_avg(endframename,endframename2)[1] - alphaterm2_avg(startframename, startframename2)[1]
    Alpha_term2 /= dtime

    #make Keplerian AM
    vkep = np.sqrt(GM1/rad_avg)
    AMkep = rad_avg*vkep*Dens

    ax0.plot(rad_avg, Dens,  lw=1.0, c=color, ls='--', label=label)
    ax0.set_ylabel(r"$\Sigma$")
    ax0.set_xlabel(r"R")
    ax0.set_xlim(0.04, np.max(rad_avg))
    ax0.set_xscale("log")
    ax0.set_ylim(0.001,100.0)
    ax0.set_yscale("log")
    ax0.legend()

    ax1.plot(rad_avg, AM, lw=1.0, c=color, ls='--')
    ax1.plot(rad_avg, AMkep, lw=1.0, c=color, ls=':', label=r'$\rm AM_{kep}$')
    ax1.set_ylabel(r"AM")
    ax1.set_xlabel(r"R")
    ax1.set_xlim(0.04, np.max(rad_avg))
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax2.plot(rad_avg, MDot, lw=1.0, c=color, ls='--')
    ax2.set_ylabel(r"$\dot{M}$")
    ax2.set_xlabel(r"R")
    ax2.set_xlim(0.04, np.max(rad_avg))
    ax2.set_xscale("log")
    ax2.set_ylim(0,0.006)

    ax3.plot(rad_avg, 2*PI*gamma*MDot/Alpha_term2, lw=1.0, c=color, ls='--')
    ax3.set_ylabel(r"$\alpha_{\rm eff}$")
    ax3.set_xlabel(r"R")
    ax3.set_xlim(0.04, np.max(rad_avg))
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylim(0.008,1.0)


#Plotting
fig = plt.figure(figsize=(12,8.5))
fig.suptitle(r"t=75-100", fontsize=15)
spec = gridspec.GridSpec(ncols=2,nrows=2)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])
ax2 = fig.add_subplot(spec[1,0])
ax3 = fig.add_subplot(spec[1,1])

plotvalues(1.1,75,100,'/LyraShared/xh2pm/BinaryDisk_prep/inflow/highres_adiabatic/stampede/gamma11_hlle_new/', 'b', r'$\gamma=1.1$', ax0, ax1, ax2, ax3)

plt.savefig("g11avg_75100.png")
#plt.show()

