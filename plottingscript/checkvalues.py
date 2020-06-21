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
gamma_gas = 1.1


#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 192
thmin = 0
thmax = 352

#load a frame that has volume
volframe = athdf('vol.athdf')

#read in pmtrack
#secondarytrack = pd.read_csv("pm_trackfile.dat",delim_whitespace=True) 

#Plotting
fig = plt.figure(figsize=(18,10.0))
spec = gridspec.GridSpec(ncols=2,nrows=2)
ax0 = fig.add_subplot(spec[0,0])
ax1 = fig.add_subplot(spec[0,1])
ax2 = fig.add_subplot(spec[1,0])
ax3 = fig.add_subplot(spec[1,1])


def deltavphi(i):
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
    deltavphi = np.array([])
    
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        v_kep = np.sqrt(GM1/rad)
        deltav = frame['vel2'][0,:,radius]-v_kep + Omega*rad
        deltavphi = np.append(deltavphi, np.sum(deltav))

                           
    return [Radius, deltavphi]

def density_avg(i):
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
    dens_avg = np.array([])
    
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        vol = volframe['user_out_var13'][0,:,radius]
        dphi = np.gradient(frame['x2v'])
        sigma =  np.sum(rho*dphi)/(2*PI)
        dens_avg = np.append(dens_avg, sigma)
                           
    return [Radius, dens_avg]

def am_avg(i):
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
    am_avg = np.array([])
    am_kep = np.array([])
    
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        sigma = rho
        vphi = frame['vel2'][0,:,radius] + Omega*rad
        vol = volframe['user_out_var13'][0,:,radius]
        dphi = np.gradient(frame['x2v'])
        am = np.sum(rad*rho*vphi*dphi)/(2*PI)
        am_avg = np.append(am_avg, am)
        am_k = np.sum(np.sqrt(GM1*rad)*sigma*dphi)/(2*PI)
        am_kep = np.append(am_kep, am_k)
        
    
    #print(np.sum(frame2['user_out_var13'][0,:,:]))
                           
    return [Radius, am_avg, am_kep]

def mdot_avg(i):
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
    mdot_avg = np.array([])

    #rvk = np.sqrt(GM1*Radius) 
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        vr = frame['vel1'][0,:,radius]
        vphi = frame['vel2'][0,:,radius] + Omega*rad
        vol = volframe['user_out_var13'][0,:,radius]
        dphi = np.gradient(frame['x2v'])
        sigma = rho
        mdot = np.sum(-rad*sigma*vr*dphi)
        mdot_avg = np.append(mdot_avg, mdot)

    return [Radius, mdot_avg]



def alphaeff_avg(i):
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
    alpha_avg = np.array([])

    #rvk = np.sqrt(GM1*Radius) 
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        press = frame['press'][0,:,radius]
        vr = frame['vel1'][0,:,radius]
        vphi = frame['vel2'][0,:,radius] + Omega*rad
        vol = volframe['user_out_var13'][0,:,radius]
        sigma = rho
        dphi = np.gradient(frame['x2v'])
        mdot = -rad*sigma*vr*dphi
        cs2 = gamma_gas*press/(rho)
        mdot_avg = np.sum(mdot)
        omega_local = vphi/rad;
        p_avg = np.sum(dphi*3*PI*sigma*cs2/omega_local)

        alpha_avg = np.append(alpha_avg, mdot_avg/p_avg)

    return [Radius, alpha_avg]


rad_avg = np.zeros_like(density_avg(150)[0]) 
Delta_vPhi = np.zeros_like(density_avg(150)[0])
Dens = np.zeros_like(density_avg(150)[0])
AM = np.zeros_like(density_avg(150)[0])
AMkep =  np.zeros_like(density_avg(150)[0])
MDot = np.zeros_like(density_avg(150)[0])
Alpha = np.zeros_like(density_avg(150)[0])

count = 0
for i in range(75,100):
    print("reading frame "+str(i))
    count += 1
    rad_avg += density_avg(i)[0]
    Delta_vPhi += deltavphi(i)[1]
    Dens += density_avg(i)[1]
    AM += am_avg(i)[1]
    AMkep += am_avg(i)[2]
    MDot += mdot_avg(i)[1]
    Alpha += alphaeff_avg(i)[1]

rad_avg /= count
Delta_vPhi /=  count
Dens /= count
AM /= count
AMkep /= count
MDot /= count
Alpha /= count


file = open('gamma11.dat','w')
file.write('rad,sigma,am,amkep,mdot,alpha\n')

for i in range(0, len(rad_avg)):
    file.write(str(rad_avg[i])+','+str(Dens[i])+','+str(AM[i])+','+str(AMkep[i])+','+str(MDot[i])+','+str(Alpha[i])+'\n')

