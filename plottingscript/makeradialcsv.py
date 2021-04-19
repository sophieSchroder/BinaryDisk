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
import pandas as pd

rc('text',usetex=False)
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

#read in pmtrack
#secondarytrack = pd.read_csv("pm_trackfile.dat",delim_whitespace=True) 

def density_avg(GM1, name, name2):

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

def am_avg(GM1, name, name2):

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

def mdot_avg(GM1, name, name2):

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


def alphaterm2_avg(GM1, name, name2):

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


def alphaeff_avg(GM1, name, name2):

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


def makecsv(GM1, gamma, start, end, dir, csvname):

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

    rad_avg = density_avg(GM1, startframename, startframename2)[0]
    Delta_vPhi = 1
    Dens = density_avg(GM1, endframename,endframename2)[1] - density_avg(GM1, startframename, startframename2)[1]
    Dens /= dtime
    AM = am_avg(GM1, endframename,endframename2)[1] - am_avg(GM1, startframename, startframename2)[1]
    AM /= dtime
    MDot = mdot_avg(GM1, endframename,endframename2)[1] - mdot_avg(GM1, startframename, startframename2)[1]
    MDot /= dtime
    Alpha = alphaeff_avg(GM1, endframename,endframename2)[1] - alphaeff_avg(GM1, startframename, startframename2)[1]
    Alpha /= dtime
    Alpha_term2 = alphaterm2_avg(GM1, endframename,endframename2)[1] - alphaterm2_avg(GM1, startframename, startframename2)[1]
    Alpha_term2 /= dtime

    #make Keplerian AM
    vkep = np.sqrt(GM1/rad_avg)
    AMkep = rad_avg*vkep*Dens

    L1 = np.max(rad_avg)
    radius = rad_avg/L1
    print("L1:", L1)

    df_data = np.array([rad_avg, radius,
                        Dens, AM, MDot, 
                        #2*PI*gamma*MDot/Alpha_term2])
                        2*PI*MDot/Alpha_term2])

    df_frame = pd.DataFrame({"rad_avg":df_data[0,:],
                             "radius":df_data[1,:],
                             "surfacedens":df_data[2,:],
                             "AM":df_data[3,:],
                             "MDot":df_data[4,:],
                             "alphaeff":df_data[5,:]})

    df_frame.to_csv(csvname, index=False)


#makecsv(0.23076923, 1.1,150,125,'./q33/g11ecc02fix/', 'q33ecc02.csv')
#makecsv(0.76923077, 1.1,150,125,'./g11ecc_hr/', 'q03ecc02.csv')
#makecsv(0.23076923, 1.1,150,125,'./q33/g11ecc035fix/', 'q33ecc035.csv')
makecsv(, 1.1, 20, 45,'./', '')


