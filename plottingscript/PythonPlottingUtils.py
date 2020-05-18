import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import ImageGrid


def read_trackfile(m1,m2,fn):
    orb=ascii.read(fn)
    
    orb['sep'] = np.sqrt(orb['x']**2 + orb['y']**2 + orb['z']**2)
    
    orb['r'] = np.array([orb['x'],orb['y'],orb['z']]).T
    orb['rhat'] = np.array([orb['x']/orb['sep'],orb['y']/orb['sep'],orb['z']/orb['sep']]).T
    
    print np.shape(orb['r'])
    orb['xcom'] = m2/m1*orb['r'][:,0]
    orb['ycom'] = m2/m1*orb['r'][:,1]
    orb['zcom'] = m2/m1*orb['r'][:,2]
    
    orb['v'] = np.array([orb['vx'],orb['vy'],orb['vz']]).T
    orb['vmag'] = np.linalg.norm(orb['v'],axis=1)
    
    if np.any(orb['vmag']!=0.0):
        orb['vhat'] = np.array([orb['vx']/orb['vmag'],orb['vy']/orb['vmag'],orb['vz']/orb['vmag']]).T
    
    
    orb['vxcom'] = m2/m1*orb['v'][:,0]
    orb['vycom'] = m2/m1*orb['v'][:,1]
    orb['vzcom'] = m2/m1*orb['v'][:,2]
    
    orb['agas1'] = np.array([orb['agas1x'],orb['agas1y'],orb['agas1z']]).T
    orb['agas2'] = np.array([orb['agas2x'],orb['agas2y'],orb['agas2z']]).T
    
    orb['rcom'] = np.array([orb['xcom'],orb['ycom'],orb['zcom']]).T
    orb['vcom'] = np.array([orb['vxcom'],orb['vy'],orb['vzcom']]).T
    
    
    return orb


def read_data(fn,orb,m1,m2,dt=1,G=1,rsoft2=0.1,level=0,
             get_cartesian=True,get_torque=False,get_energy=False,
             x1_min=None,x1_max=None,
             x2_min=None,x2_max=None,
             x3_min=None,x3_max=None,
             profile_file="hse_profile.dat"):
    """ Read spherical data and reconstruct cartesian mesh for analysis/plotting """
  
    
    print "read_data...reading file",fn
    
    
    d = ar.athdf(fn,level=level,subsample=True,
                 x1_min=x1_min,x1_max=x1_max,
                 x2_min=x2_min,x2_max=x2_max,
                 x3_min=x3_min,x3_max=x3_max) # approximate arrays by subsampling if level < max
    print " ...file read, constructing arrays"
    
    # current time
    t = d['Time']
    # get properties of orbit
    rcom,vcom = rcom_vcom(orb,t)
    x2,y2,z2 = pos_secondary(orb,t)
    
    
    # MAKE grid based coordinates
    d['gx1v'] = np.zeros_like(d['rho'])
    for i in range((d['rho'].shape)[2]):
        d['gx1v'][:,:,i] = d['x1v'][i]
    
    d['gx2v'] = np.zeros_like(d['rho'])
    for j in range((d['rho'].shape)[1]):
        d['gx2v'][:,j,:] = d['x2v'][j]

    d['gx3v'] = np.zeros_like(d['rho'])
    for k in range((d['rho'].shape)[0]):
        d['gx3v'][k,:,:] = d['x3v'][k]
    
    
    ####
    # GET THE VOLUME 
    ####
    
    ## dr, dth, dph
    d1 = d['x1f'][1:] - d['x1f'][:-1]
    d2 = d['x2f'][1:] - d['x2f'][:-1]
    d3 = d['x3f'][1:] - d['x3f'][:-1]
    
    # grid based versions
    gd1 = np.zeros_like(d['rho'])
    for i in range((d['rho'].shape)[2]):
        gd1[:,:,i] = d1[i]
    
    gd2 = np.zeros_like(d['rho'])
    for j in range((d['rho'].shape)[1]):
        gd2[:,j,:] = d2[j]

    gd3 = np.zeros_like(d['rho'])
    for k in range((d['rho'].shape)[0]):
        gd3[k,:,:] = d3[k]
    
    # AREA / VOLUME 
    sin_th = np.sin(d['gx2v'])
    d['dA'] = d['gx1v']**2 * sin_th * gd2*gd3
    d['dvol'] = d['dA'] * gd1
    
    # free up d1,d2,d3
    del d1,d2,d3
    del gd1,gd2,gd3
    
    
    ### 
    # CARTESIAN VALUES
    ###
    if(get_cartesian or get_torque or get_energy):
        # angles
        cos_th = np.cos(d['gx2v'])
        sin_ph = np.sin(d['gx3v'])
        cos_ph = np.cos(d['gx3v']) 
        
        # cartesian coordinates
        d['x'] = d['gx1v'] * sin_th * cos_ph 
        d['y'] = d['gx1v'] * sin_th * sin_ph 
        d['z'] = d['gx1v'] * cos_th
    
        # cartesian velocities
        d['vx'] = sin_th*cos_ph*d['vel1'] + cos_th*cos_ph*d['vel2'] - sin_ph*d['vel3'] 
        d['vy'] = sin_th*sin_ph*d['vel1'] + cos_th*sin_ph*d['vel2'] + cos_ph*d['vel3'] 
        d['vz'] = cos_th*d['vel1'] - sin_th*d['vel2']  

        del cos_th, sin_th, cos_ph, sin_ph
    
    
    if(get_torque):
        
        # define grav forces
        dist2 = np.sqrt( (d['x']-x2)**2 + (d['y']-y2)**2 + (d['z']-z2)**2 )
        dist1c = d['gx1v']**3
    
        soft_grav = fspline(dist2,rsoft2)
        
        fdens2x = G*m2*d['rho']*soft_grav * (d['x']-x2)
        fdens2y = G*m2*d['rho']*soft_grav * (d['y']-y2)
        fdens2z = G*m2*d['rho']*soft_grav * (d['z']-z2)

        fdens1x = G*m1*d['rho']/dist1c * d['x']
        fdens1y = G*m1*d['rho']/dist1c * d['y']
        fdens1z = G*m1*d['rho']/dist1c * d['z']

        del dist1c,dist2
        
        d['torque_dens_2_z'] = (d['x']-rcom[0])*fdens2y - (d['y']-rcom[1])*fdens2x
        d['torque_dens_1_z'] = (d['x']-rcom[0])*fdens1y - (d['y']-rcom[1])*fdens1x
        
        del fdens2x,fdens2y,fdens2z
        del fdens1x,fdens1y,fdens1z

    if(get_energy):
       #energy
        dist2 = np.sqrt( (d['x']-x2)**2 + (d['y']-y2)**2 + (d['z']-z2)**2 )
        # should update with the real potential 
        epot = -G*m1*d['rho']/d['gx1v'] - G*m2*d['rho']*pspline(dist2,rsoft2)
        ek = 0.5*d['rho']*((d['vx']-vcom[0])**2 +
                           (d['vy']-vcom[1])**2 + 
                           (d['vz']-vcom[2])**2)
        ei = d['press']/(5/3.-1)
        d['etot'] = epot + ei + ek
        
        del dist2,epot,ek,ei
    
    return d


def rcom_vcom(orb,t):
    """pass a pm_trackfile.dat that has been read, time t"""
    
    #if(corotating):
        
    #### pm_trackfile has no rcom and vcom 
    rcom =  np.array([np.interp(t,orb['time'],orb['rcom'][:,0]),
                  np.interp(t,orb['time'],orb['rcom'][:,1]),
                  np.interp(t,orb['time'],orb['rcom'][:,2])])
    vcom =  np.array([np.interp(t,orb['time'],orb['vcom'][:,0]),
                  np.interp(t,orb['time'],orb['vcom'][:,1]),
                  np.interp(t,orb['time'],orb['vcom'][:,2])])

    return rcom,vcom

def pos_secondary(orb,t):
    x2 = np.interp(t,orb['time'],orb['x'])
    y2 = np.interp(t,orb['time'],orb['y'])
    z2 = np.interp(t,orb['time'],orb['z'])
    return x2,y2,z2

    
# individual grav force 
def fspline(r,eps):
    """Hernquist & Katz 1989, Eq A2 """
    u=r/eps
    condlist = [u<1,u<2,u>=2]
    resultlist = [1/eps**3 * (4./3. - 1.2*u**2 + 0.5*u**3),
                  1/r**3 * (-1./15. + 8./3.*u**3 - 3.*u**4 + 1.2*u**5 - 1./6.*u**6),
                  1/r**3]
    return np.select(condlist,resultlist)


def pspline(r,eps):
    """Hernquist & Katz 1989, Eq A1 """
    u=r/eps
    condlist = [u<1,u<2,u>=2]
    resultlist = [-2./eps * (1./3.*u**2 -0.15*u**4 + 0.05*u**5) + 1.2*eps,
                  1./(15.*r) - 1/eps*( 4./3.*u**2 - u**3 + 0.3*u**4 -1./30.*u**5) + 1.6*eps,
                  1./r]
    return np.select(condlist,resultlist)



def get_plot_array_midplane(arr):
    return np.append(arr,[arr[0]],axis=0)
    
def get_plot_array_edge(arr):
    print np.shape(arr)
    
    print np.shape([arr[0]])
    print np.shape(np.array([arr[:,0]]).T)
    
    
    step1 = np.append(arr,np.array([arr[:,0]]).T,axis=1)
    print np.shape(step1)
    step2 = np.zeros(shape=(1,1))
    step3 = np.append(step1,np.append(step2,[arr[0]],axis=1),axis=0)
    
    return step3


