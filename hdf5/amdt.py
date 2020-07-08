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
rmax = 384
thmin = 0
thmax = 704


#read in pmtrack
#secondarytrack = pd.read_csv("pm_trackfile.dat",delim_whitespace=True) 

#Plotting
fig = plt.figure(figsize=(10,7.5))
spec = gridspec.GridSpec(ncols=1,nrows=1)
ax0 = fig.add_subplot(spec[:])

jobid = "thindisk"

def readam(i):
    #change file names here accordingly
    if i < 10:
        name = jobid+'.out1.0000'+str(i)+'.athdf'
        name2 = jobid+'.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = jobid+'.out1.00'+str(i)+'.athdf'
        name2 = jobid+'.out2.00'+str(i)+'.athdf'
    else:
        name = jobid+'.out1.000'+str(i)+'.athdf'
        name2 = jobid+'.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    AM = np.array([])

    for radius in range(rmin,rmax):
        
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        v_kep = np.sqrt(GM1/rad)
        deltav = frame['vel2'][0,:,radius]-v_kep +Omega*rad # or not adding Omega*rad?
        #am = rho*rad*deltav
        #am_az = np.average(am)*rad
        am = np.sum(frame2['user_out_var9'][0,:,radius])
        vol = np.sum(frame2['user_out_var13'][0,:,radius])
        am_az = am*rad #/vol
        AM = np.append(AM, am_az)
        
    return [Radius,AM,time]


def damdt(index): #backward gradient
    
    data_now = readam(index)
    data_before = readam(index-1)
    
    rad_now = data_now[0]
    rad_before = data_before[0]

    AM_now = data_now[1]
    AM_before = data_before[1]

    t_now = data_now[2]
    t_before = data_before[2]

    dAM = AM_now - AM_before
    dt = t_now - t_before
    dAMdt = dAM/dt
    

    return [rad_now, dAMdt]


def ammdot(i):
    #change file names here accordingly
    if i < 10:
        name = jobid+'.out1.0000'+str(i)+'.athdf'
        name2 = jobid+'.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = jobid+'.out1.00'+str(i)+'.athdf'
        name2 = jobid+'.out2.00'+str(i)+'.athdf'
    else:
        name = jobid+'.out1.000'+str(i)+'.athdf'
        name2 = jobid+'.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)

    time = frame['Time']
    
    Radius = frame['x1v']
    AMMdot = np.array([])

    rvk =  np.sqrt(GM1*frame['x1v'])
    drvkdr = np.gradient(rvk)/np.gradient(Radius)

    for radius in range(rmin,rmax):
        
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        v_kep = np.sqrt(GM1/rad)
        deltav = frame['vel2'][0,:,radius]-v_kep +Omega*rad
        vr = frame['vel1'][0,:,radius]
        
        #mdot_arr = -rho*rad*vr
        #mdot = np.average(mdot_arr)

        mdot = np.sum(frame2['user_out_var10'][0,:,radius])
        AMMdot = np.append(AMMdot, mdot)

    AMMdot *= drvkdr
        
    return [Radius,AMMdot]



def amfh(i):
    if i < 10:
        name = jobid+'.out1.0000'+str(i)+'.athdf'
        name2 = jobid+'.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = jobid+'.out1.00'+str(i)+'.athdf'
        name2 = jobid+'.out2.00'+str(i)+'.athdf'
    else:
        name = jobid+'.out1.000'+str(i)+'.athdf'
        name2 = jobid+'.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)
    time = frame['Time']
    
    Radius = frame['x1v']
    amfh = np.array([])
    
    for radius in range(rmin,rmax):
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        v_kep = np.sqrt(GM1/rad)
        deltav = frame['vel2'][0,:,radius]-v_kep + Omega*rad
        vr = frame['vel1'][0,:,radius]        

        #FH = rho*vr*deltav
        #FH = np.average(FH)
        FH = np.sum(frame2['user_out_var11'][0,:,radius])
        amfh = np.append(amfh, FH*rad*rad)

    AMFH = -np.gradient(amfh)/np.gradient(Radius)
                           
    return [Radius, AMFH]

    
def torque(i):
    #change file names here accordingly
    if i < 10:
        name = jobid+'.out1.0000'+str(i)+'.athdf'
        name2 = jobid+'.out2.0000'+str(i)+'.athdf'
    elif i >=100:
        name = jobid+'.out1.00'+str(i)+'.athdf'
        name2 = jobid+'.out2.00'+str(i)+'.athdf'
    else:
        name = jobid+'.out1.000'+str(i)+'.athdf'
        name2 = jobid+'.out2.000'+str(i)+'.athdf'

    frame = athdf(name)
    frame2 = athdf(name2)

    time = frame['Time']
    
    Radius = frame['x1v']
    amfh = np.array([])

    trr = np.array([])
    trr2 = np.array([])

    #pmtrack_index = i-1
    #x2 = secondarytrack.x[pmtrack_index]
    #y2 = secondarytrack.y[pmtrack_index]
    #z2 = secondarytrack.z[pmtrack_index]
    
    for radius in range(rmin,rmax):
        fext = frame2['user_out_var7'][0,:,radius]
        tr = Radius[radius]*fext
        #torque = np.average(tr)
        torque = np.sum(frame2['user_out_var12'][0,:,radius])
        trr = np.append(trr, torque)
        '''
        torq_arr = np.array([])
        
        for th in range(thmin, thmax):
        
            rad = Radius[radius]
            theta = frame['x2v'][th]
            #z_cyl = frame['x3v'][radius]

            x = rad*np.cos(theta)
            y = rad*np.sin(theta)
            #z = z_cyl

            r1 = np.sqrt(x*x+y*y) #+z*z)
            r12 = np.sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y))
            f_m2_x = (x2-x)*GM2*frame['rho'][0,th,radius]/pow(r12,3)
            f_m2_y = (y2-y)*GM2*frame['rho'][0,th,radius]/pow(r12,3)
            torq = x*f_m2_y - y*f_m2_x

            torq_arr = np.append(torq_arr, torq)
        
        trr2 = np.append(trr2, np.average(torq_arr))
        '''
    Tr = trr*Radius
    #Tr2 = trr2*Radius
    return [Radius, Tr]#,  Tr2]

'''
dAMdt = damdt(400)
#AM_Mdot = ammdot(154)
#AM_FH = amfh(154)
#Tr = torque(154)


ax0.plot(dAMdt[0], dAMdt[1], label=r"$AM_{t}(R)$")
#ax0.plot(AM_Mdot[0], AM_Mdot[1], label=r"$AM_{\dot{M}}(R)$")
#ax0.plot(AM_FH[0], AM_FH[1], label=r"$AM_{FH}(R)$")
#ax0.plot(Tr[0], Tr[1], label=r"$T(R)$")
#ax0.plot(Tr[0], Tr[2], label="calculated T(R)")
#ax0.plot(AM_Mdot[0], AM_Mdot[1]+AM_FH[1]+Tr[1], ls='--', label=r"$AM_{\dot{M}}(R)+AM_{FH}(R)+T(R)$")
#print(np.shape(AM_Mdot[0]), np.shape(AM_FH[0]), np.shape(Tr[0]))

#ax0.set_xlim(np.min(Tr[0]), np.max(Tr[0]))
ax0.set_xlabel(r"R")

plt.show()

'''


rad_avg = np.zeros_like(damdt(400)[0]) #np.array([])
dAMdt_avg = np.zeros_like(damdt(400)[1])
AM_Mdot_avg = np.zeros_like(damdt(400)[1])
AM_FH_avg = np.zeros_like(damdt(400)[1])
Tr_avg = np.zeros_like(damdt(400)[1])

count = 0
for i in range(395,400):
    count += 1
    rad_avg += damdt(i)[0]
    dAMdt_avg += damdt(i)[1]
    AM_Mdot_avg += ammdot(i)[1]
    AM_FH_avg += amfh(i)[1]
    Tr_avg += torque(i)[1]

rad_avg /= count #np.append(rad_avg, dAMdt[0])
dAMdt_avg /= count #np.append(dAMdt_avg, dAMdt[1])
AM_Mdot_avg /= count #np.append(AM_Mdot_avg, AM_Mdot[1])
AM_FH_avg /= count #np.append(AM_Mdot_avg, AM_FH[1])
Tr_avg /= count #np.apppend(Tr_avg, Tr[1])

ax0.plot(rad_avg, dAMdt_avg, label=r"$AM_{t}(R)$", lw=1.0, c='b')
ax0.plot(rad_avg, AM_Mdot_avg, lw=1.0, label=r"$AM_{\dot{M}}(R)$")
ax0.plot(rad_avg, AM_FH_avg, lw=1.0, label=r"$AM_{FH}(R)$")
ax0.plot(rad_avg, Tr_avg, lw=1.0, label=r"$T(R)$")
ax0.plot(rad_avg, dAMdt_avg-AM_Mdot_avg-Tr_avg, ls='--', lw=1.0, c='b')

#ax0.set_xlim(np.min(rad_avg[0]), np.max(rad_avg[0]))
ax0.set_xlabel(r"R")
#ax0.set_ylim(-1.0e-5, 1.0e-5)


plt.legend()
plt.show()
#plt.savefig("AMdt.png")

