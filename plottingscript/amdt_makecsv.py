#!/usr/bin/env python
import numpy as np
from athena_read import athdf #need athena_read.py, it's in ~/code/vis/python
import numpy as np
import matplotlib.pyplot as plt
from athena_read import athdf
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


#reading angular momentum
def readam(frame2):

    time = frame2['Time']
    
    Radius = frame2['x1v']

    AM = np.sum(frame2['user_out_var6'][0,:,:], axis=0)
    AM_net = np.sum(frame2['user_out_var11'][0,:,:], axis=0)
    AM_nc = np.sum(frame2['user_out_var34'][0,:,:], axis=0)
        
    return [Radius, AM, time, AM_net, AM_nc]

def intammdot(frame2):

    time = frame2['Time']
    
    Radius = frame2['x1v']
    dR = np.gradient(Radius)
    AMMdot = np.sum(frame2['user_out_var16'][0,:,:], axis=0)

    return [Radius,AMMdot]


def intamfh(frame2):

    time = frame2['Time']
    
    Radius = frame2['x1v']
    dR = np.gradient(Radius)
    amfh = np.sum(frame2['user_out_var17'][0,:,:], axis=0)
    amth_net = np.sum(frame2['user_out_var18'][0,:,:], axis=0)
    amth_nc = np.sum(frame2['user_out_var37'][0,:,:], axis=0)

    return [Radius, amfh, amth_net, amth_nc]



def inttorque(frame2):

    time = frame2['Time']
    
    Radius = frame2['x1v']
    dR = np.gradient(Radius)
    trr = np.sum(frame2['user_out_var19'][0,:,:], axis=0)
    trr2 = np.sum(frame2['user_out_var20'][0,:,:], axis=0)
        
    return [Radius, trr,  trr2]


def makecsv(start_idx, end_idx, csvname):
    if start_idx< 10:
        name_s = 'BDstream.out2.0000'+str(start_idx)+'.athdf'
    elif start_idx >=100:
        name_s = 'BDstream.out2.00'+str(start_idx)+'.athdf'
    else:
        name_s = 'BDstream.out2.000'+str(start_idx)+'.athdf'

    if end_idx < 10:
        name_e = 'BDstream.out2.0000'+str(end_idx)+'.athdf'
    elif end_idx >=100:
        name_e = 'BDstream.out2.00'+str(end_idx)+'.athdf'
    else:
        name_e = 'BDstream.out2.000'+str(end_idx)+'.athdf'

    start = athdf(name_s)
    end = athdf(name_e)

    rad_avg = readam(end)[0]
    dAM_r = readam(end)[1]-readam(start)[1]
    dAM_nr = readam(end)[4]-readam(start)[4]
    AMMdot = intammdot(end)[1] - intammdot(start)[1]
    AMTH_r = intamfh(end)[1] - intamfh(start)[1]
    AMTH_nr = intamfh(end)[3] - intamfh(start)[3]
    TR_nr = inttorque(end)[2] - inttorque(start)[2]
    TR_r  = inttorque(end)[1] - inttorque(start)[1]
    TRdiff = TR_r - TR_nr
    dtime = readam(end)[2]-readam(start)[2]

    dAM_r /= dtime
    dAM_nr /= dtime
    AMMdot /= dtime
    AMTH_r /= dtime
    AMTH_nr /= dtime
    TR_r /= dtime
    TR_nr /= dtime
    TRdiff /= dtime

    df_data = np.array([rad_avg, dAM_r, dAM_nr, AMMdot,\
                        AMTH_r, AMTH_nr, TR_nr, TR_r, \
                        TRdiff, dtime], dtype=object)

    df_frame = pd.DataFrame({"rad_avg":df_data[0],\
                             "dAM_r":df_data[1],\
                             "dAM_nr":df_data[2],\
                             "AMMdot":df_data[3],\
                             "AMTH_r":df_data[4],\
                             "AMTH_nr":df_data[5],\
                             "TR_nr":df_data[6],\
                             "TR_r":df_data[7],\
                             "TRdiff":df_data[8]})

    df_frame.to_csv(csvname, index=False)

makecsv(100,75,"diskscale01_amdt.csv")

