#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from athena_read import athdf
import numpy as np
from time import time
import multiprocessing as mp
import matplotlib.pyplot as plt
from athena_read import hst
import matplotlib.gridspec as gridspec
from matplotlib import animation
from matplotlib import rc 

rc('text',usetex=True)
rc('font',family='serif',size=12)



#reading history file
hst = hst('BDstream.hst')

#Plotting
fig = plt.figure(figsize=(10,7.5))
spec = gridspec.GridSpec(ncols=1,nrows=2)
ax0 = fig.add_subplot(spec[0,:])
ax1 = fig.add_subplot(spec[1,:])

ax0.plot(hst['time'][1:], -hst['massfluxix1'][1:], label=r'$\dot{M}_{\rm inner}$')
ax0.plot(hst['time'][1:], -hst['massfluxox1'][1:], label=r'$\dot{M}_{\rm outer}$')
#ax0.set_yscale('log')
ax0.legend(loc='upper right')
ax0.set_ylabel(r'$\dot{M}$')

massinflow = hst['massfluxix1']-hst['massfluxox1']
cumsum = np.cumsum(massinflow)
dt = 0.1
dv = (0.6/128)*(2*np.pi/256)
print(dv)

ax1.plot(hst['time'][1:], hst['mass'][1:])
#ax1.plot(hst['time'][1:], cumsum[1:]*dv)
ax1.set_xlabel(r'$t/t_{0}$')
ax1.set_ylabel(r'$M_{\rm tot}/\Sigma_{0}A_{0}$')

plt.show()
