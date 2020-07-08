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
gamma_gas = 1.3


#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 192
thmin = 0
thmax = 352

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
        mdot =  frame2['user_out_var23'][0,:,radius]
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
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        press = frame['press'][0,:,radius]
        vr = frame['vel1'][0,:,radius]
        vphi = frame['vel2'][0,:,radius] + Omega*rad
        sigma = rho
        dphi = np.gradient(frame['x2v'])
        mdot = -rad*sigma*vr*dphi
        cs2 = gamma_gas*press/(rho)
        mdot_avg = np.sum(mdot)
        omega_local = vphi/rad;
        p_avg = np.sum(dphi*3*PI*sigma*cs2/omega_local)

        alpha_eff = frame2['user_out_var22'][0,:,radius]
        alpha_eff  = np.sum(alpha_eff)/(2*PI)
        alpha_avg = np.append(alpha_avg, alpha_eff)

    return [Radius, alpha_avg]



startframename = 'BDstream.out1.00050.athdf'
endframename = 'BDstream.out1.00100.athdf'
startframename2 = 'BDstream.out2.00050.athdf'
endframename2 = 'BDstream.out2.00100.athdf'
starttime = athdf(startframename)['Time']
endtime = athdf(endframename)['Time']
dtime = endtime-starttime

rad_avg = density_avg(startframename, startframename2)[0]
Delta_vPhi = 1
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

startframename_ = './lowfloor/BDstream.out1.00050.athdf'
endframename_ = './lowfloor/BDstream.out1.00100.athdf'
startframename2_ = './lowfloor/BDstream.out2.00050.athdf'
endframename2_ = './lowfloor/BDstream.out2.00100.athdf'
starttime_ = athdf(startframename_)['Time']
endtime_ = athdf(endframename_)['Time']
dtime_ = endtime-starttime

rad_avg_ = density_avg(startframename_, startframename2_)[0]
Dens_ = density_avg(endframename_,endframename2_)[1] - density_avg(startframename_, startframename2_)[1]
Dens_ /= dtime_
AM_ = am_avg(endframename_,endframename2_)[1] - am_avg(startframename_, startframename2_)[1]
AM_ /= dtime_
MDot_ = mdot_avg(endframename_,endframename2_)[1] - mdot_avg(startframename_, startframename2_)[1]
MDot_ /= dtime_
Alpha_ = alphaeff_avg(endframename_,endframename2_)[1] - alphaeff_avg(startframename_, startframename2_)[1]
Alpha_term2_ = alphaterm2_avg(endframename_,endframename2_)[1] - alphaterm2_avg(startframename_, startframename2_)[1]
Alpha_term2_ /= dtime_

#Plotting
fig = plt.figure(figsize=(12,8.5))
fig.suptitle(r"$\gamma=1.1, t=50-100$", fontsize=15)
spec = gridspec.GridSpec(ncols=2,nrows=2)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])
ax2 = fig.add_subplot(spec[1,0])
ax3 = fig.add_subplot(spec[1,1])

ax0.plot(rad_avg, Dens,  lw=1.0, c='b', ls='--', label=r"original floor")
ax0.plot(rad_avg_, Dens_,  lw=1.0, c='b', ls=':', label=r'lower floor')
ax0.set_ylabel(r"$\Sigma$")
ax0.set_xlabel(r"R")
ax0.set_xlim(0.04, np.max(rad_avg))
ax0.set_xscale("log")
ax0.set_ylim(0.001,100.0)
ax0.set_yscale("log")
ax0.legend()

ax1.plot(rad_avg, AM, lw=1.0, c='b', ls='--')
ax1.plot(rad_avg_, AM_, lw=1.0, c='b', ls=':')
ax1.set_ylabel(r"AM")
ax1.set_xlabel(r"R")
ax1.set_xlim(0.04, np.max(rad_avg))
ax1.set_xscale("log")
ax1.set_yscale("log")

ax2.plot(rad_avg, MDot, lw=1.0, c='b', ls='--')
ax2.plot(rad_avg_, MDot_, lw=1.0, c='b', ls=':')
ax2.set_ylabel(r"$\dot{M}$")
ax2.set_xlabel(r"R")
ax2.set_xlim(0.04, np.max(rad_avg))
ax2.set_xscale("log")

ax3.plot(rad_avg, MDot/Alpha_term2, lw=1.0, c='b', ls='--')
ax3.plot(rad_avg_, MDot_/Alpha_term2_, lw=1.0, c='b', ls=':')
ax3.set_ylabel(r"$\alpha_{\rm eff}$")
ax3.set_xlabel(r"R")
ax3.set_xlim(0.04, np.max(rad_avg))
ax3.set_xscale("log")
ax3.set_yscale("log")
#ax3.set_ylim(0.001,1.0)

plt.savefig("gamma11_avg_floor.png")
#plt.show()
