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

rc('text',usetex=True)
rc('font',family='serif',size=12)

GM1 = 0.7692307692307692 

#the range in R direction we want to plot, it's the cell number, not physical distance, change them as you want
rmin = 0
rmax = 128


#Plotting
fig = plt.figure(figsize=(10,7.5))
spec = gridspec.GridSpec(ncols=1,nrows=1)
ax0 = fig.add_subplot(spec[:])

def readam(ax, i):
    #change file names here accordingly
    if i < 10:
        name = 'BDstream.out1.0000'+str(i)+'.athdf'
    elif i >=100:
        name = 'BDstream.out1.00'+str(i)+'.athdf'
    else:
        name = 'BDstream.out1.000'+str(i)+'.athdf'

    frame = athdf(name)
    time = frame['Time']
    
    Radius = frame['x1v']
    AM = np.array([])

    for radius in range(rmin,rmax):
        
        rad = frame['x1v'][radius]
        rho = frame['rho'][0,:,radius]
        v_kep = np.sqrt(GM1/rad)
        am = user_out_var9
        
        AM = np.append(AM, am_az)
        
    ax.plot(Radius, AM, label=r'time='+str(time))


readam(ax0,140)
readam(ax0,150)
readam(ax0,155)
plt.legend()
plt.show()
