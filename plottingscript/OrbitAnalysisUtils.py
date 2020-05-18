# normal stuff
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
    #orb['lgoz'] = np.cumsum( np.gradient(orb['time']) * orb['ldoz'] )
    #orb['ltz'] = orb['lpz'] + orb['lgz'] + orb['lgoz']
    
    orb['sep'] = np.sqrt(orb['x']**2 + orb['y']**2 + orb['z']**2)
    
    orb['r'] = np.array([orb['x'],orb['y'],orb['z']]).T
    orb['rhat'] = np.array([orb['x']/orb['sep'],orb['y']/orb['sep'],orb['z']/orb['sep']]).T
    
    orb['v'] = np.array([orb['vx'],orb['vy'],orb['vz']]).T
    orb['vmag'] = np.linalg.norm(orb['v'],axis=1)
    orb['vhat'] = np.array([orb['vx']/orb['vmag'],orb['vy']/orb['vmag'],orb['vz']/orb['vmag']]).T
    
    orb['agas1'] = np.array([orb['agas1x'],orb['agas1y'],orb['agas1z']]).T
    orb['agas2'] = np.array([orb['agas2x'],orb['agas2y'],orb['agas2z']]).T
    
    orb['rcom'] = np.array([orb['x'],orb['y'],orb['z']]).T * m2/(m1+m2)
    orb['vcom'] = np.array([orb['vx'],orb['vy'],orb['vz']]).T * m2/(m1+m2)
        
    # np.array([orb['xcom'],orb['ycom'],orb['zcom']]).T
    # np.array([orb['vxcom'],orb['vycom'],orb['vzcom']]).T
    
    F12 = - m1*m2/orb['sep']**2
    # aceel of 1 by 2
    orb['a21'] = np.array([-F12/m1*orb['x']/orb['sep'],
                           -F12/m1*orb['y']/orb['sep'],
                           -F12/m1*orb['z']/orb['sep']]).T
    # accel of 2 by 1
    orb['a12'] = np.array([F12/m2*orb['x']/orb['sep'],
                           F12/m2*orb['y']/orb['sep'],
                           F12/m2*orb['z']/orb['sep']]).T
    
    #mu = 0.631686 + 0.339421 + 0.3
    #orb['E'] = orb['vmag']**2 / 2. - mu/orb['sep']
    #orb['a'] = - mu / (2*orb['E'])
    
    
    return orb


def get_t1(orb,skip=1):
    sel = orb['sep']<1.5
    return np.interp(1.0,np.flipud(orb[sel]['sep'][::skip]),np.flipud(orb[sel]['time'][::skip]) )




def get_Omega_env_dist(fn,dv=0.05,G=1,rho_thresh=1.e-2,level=2):
    """ Get the mass-average Omega within r<1 """
    d = ar.athdf(fn,level=level)

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
    
    
    ## dr, dth, dph
    d['d1'] = d['x1f'][1:] - d['x1f'][:-1]
    d['d2'] = d['x2f'][1:] - d['x2f'][:-1]
    d['d3'] = d['x3f'][1:] - d['x3f'][:-1]

    # grid based versions
    d['gd1'] = np.zeros_like(d['rho'])
    for i in range((d['rho'].shape)[2]):
        d['gd1'][:,:,i] = d['d1'][i]
    
    d['gd2'] = np.zeros_like(d['rho'])
    for j in range((d['rho'].shape)[1]):
        d['gd2'][:,j,:] = d['d2'][j]

    d['gd3'] = np.zeros_like(d['rho'])
    for k in range((d['rho'].shape)[0]):
        d['gd3'][k,:,:] = d['d3'][k]
    
    
    # VOLUME 
    d['dvol'] = d['gx1v']**2 * np.sin(d['gx2v']) * d['gd1']* d['gd2']* d['gd3']
    # cell masses
    d['dm'] = d['rho']*d['dvol']


    select = (d['rho']>rho_thresh) #(d['gx1v']<1.0)
    # Omega = vphi/(r sin th)
    vpf =  (d['vel3'][select] / (d['gx1v'][select]* np.sin(d['gx2v'][select])) ).flatten()
    #vpf =  d['vel3'][select].flatten()
    dmf = d['dm'][select].flatten()
    
    mybins = np.arange(-0.1,0.1,dv)
    mydist = np.histogram(vpf,weights=dmf,bins=mybins)
    GMtot=G*np.sum(mydist[0])
    
    print "Total mass GM = ",G*np.sum(dmf)
    print "Total mass GM (distribution) = ", GMtot
    
    return mydist, np.average(vpf,weights=dmf)


def read_data(fn,orb,m1,m2,G=1,rsoft2=0.1,level=0,
             get_cartesian=True,get_torque=False,get_energy=False,
             rotate_frame=False,
             x1_min=None,x1_max=None,
             x2_min=None,x2_max=None,
             x3_min=None,x3_max=None,
             profile_file="hse_profile.dat",
             gamma=5./3.  ):
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
    # not necessary for wind setup
    #rcom = np.array([0.5,0.0,0.0]).T
    #vcom = np.array([0.0,0.0,0.0]).T
    #x2,y2,z2 = np.array([1.0,0.0,0.0]).T
    
    # MAKE grid based coordinates
    d['gx1v'] = np.zeros_like(d['rho'])
    for i in range((d['vel1'].shape)[2]):
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
    gd1 = np.zeros_like(d['vel1'])
    for i in range((d['vel1'].shape)[2]):
        gd1[:,:,i] = d1[i]
    
    gd2 = np.zeros_like(d['vel1'])
    for j in range((d['vel1'].shape)[1]):
        gd2[:,j,:] = d2[j]

    gd3 = np.zeros_like(d['vel1'])
    for k in range((d['vel1'].shape)[0]):
        gd3[k,:,:] = d3[k]
    
    # AREA / VOLUME 
    sin_th = np.sin(d['gx2v'])
    d['dA'] = d['gx1v']**2 * sin_th * gd2*gd3
    d['dvol'] = d['dA'] * gd1

    # save for timestep
    d['gd1'] = gd1
    d['gd2'] = gd2
    d['gd3'] = gd3
    
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
        
        if(rotate_frame):
            
            # cartesian coordinates
            d['x'] = d['gx1v'] * sin_th * cos_ph 
            d['y'] = d['gx1v'] * sin_th * sin_ph 
            d['z'] = d['gx1v'] * cos_th
    
            # cartesian velocities
            d['vx'] = sin_th*cos_ph*d['vel1'] + cos_th*cos_ph*d['vel2'] - sin_ph*d['vel3'] 
            d['vy'] = sin_th*sin_ph*d['vel1'] + cos_th*sin_ph*d['vel2'] + cos_ph*d['vel3'] 
            d['vz'] = cos_th*d['vel1'] - sin_th*d['vel2']  
            
            
            # rotate coordinates 
            # omega = 1.0 :  transform corotating coordinates to orbit frame
            # omega = -1.0:  transform orbit coordinates to corotating frame

            # get vphi in vx, vy
            vphi = omega_orb*d['gx1v']*sin_th

            d['vxrot'] = - sin_ph*vphi
            d['vyrot'] =  cos_ph*vphi
            # d['vzrot'] =  

            # rotate coordinates
            d['xorb'] = np.cos(omega_orb*d['Time'])*d['x'] - np.sin(omega_orb*d['Time'])*d['y']
            d['yorb'] = np.sin(omega_orb*d['Time'])*d['x'] + np.cos(omega_orb*d['Time'])*d['y']
            d['zorb'] = d['z'] 

            # rotate velocities
            d['vxorb'] = np.cos(omega_orb*d['Time'])*(d['vx']+d['vxrot']) - np.sin(omega_orb*d['Time'])*(d['vy']+d['vyrot']) 
            d['vyorb'] = np.sin(omega_orb*d['Time'])*(d['vx']+d['vxrot']) + np.cos(omega_orb*d['Time'])*(d['vy']+d['vyrot'])
            d['vzorb'] = d['vz']

        else:
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
        #fdens2z = G*m2*d['rho']*soft_grav * (d['z']-z2)

        fdens1x = G*m1*d['rho']/dist1c * d['x']
        fdens1y = G*m1*d['rho']/dist1c * d['y']
        #fdens1z = G*m1*d['rho']/dist1c * d['z']

        del dist1c,dist2

        d['torque_dens_2_z'] = (x2-rcom[0])*fdens2y - (y2-rcom[1])*fdens2x
        d['torque_dens_1_z'] = (-rcom[0])*fdens1y - (-rcom[1])*fdens1x
        
        del fdens2x,fdens2y #,fdens2z
        del fdens1x,fdens1y #,fdens1z

    if(get_energy):
        hse_prof = ascii.read(profile_file,
                              names=['r','rho','p','m'])
        #energy
        dist2 = np.sqrt( (d['x']-x2)**2 + (d['y']-y2)**2 + (d['z']-z2)**2 )
        M1r = np.interp(d['gx1v'],hse_prof['r'],hse_prof['m'])
        # should update with the real potential
        d['epotg'] = -G*(M1r-m1)*d['rho']/d['gx1v']
        d['epotp'] = -G*m1*d['rho']/d['gx1v'] - G*m2*d['rho']*pspline(dist2,rsoft2)
        d['ek'] = 0.5*d['rho']*((d['vx']-vcom[0])**2 +
                           (d['vy']-vcom[1])**2 + 
                           (d['vz']-vcom[2])**2)
        d['ei'] = d['press']/(gamma-1)
        d['etot'] = d['epotg'] + d['epotp'] + d['ei'] + d['ek']
        d['h'] = gamma*d['press']/((gamma-1)*d['rho'])
        d['bern'] = (d['etot']+d['press'])/d['rho']
        
        del hse_prof,dist2,M1r
    
    return d
    
# get time from filename
def time_fn(fn,dt=1):
    """ returns a float time from parsing the fn string (assumes dt=1 in fn convention)"""
    return dt*float(fn[-11:-6])


def rcom_vcom(orb,t):
    """pass a pm_trackfile.dat that has been read, time t"""
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
    resultlist = [-2./eps * (1./3.*u**2 -0.15*u**4 + 0.05*u**5) + 7./(5.*eps),
                  -1./(15.*r) - 1/eps*( 4./3.*u**2 - u**3 + 0.3*u**4 -1./30.*u**5) + 8./(5.*eps),
                  1./r]
    return np.select(condlist,resultlist)





def get_plot_array_midplane(arr):
    return np.append(arr,[arr[0]],axis=0)



from scipy.signal import argrelextrema

class makebinary:
    """assumes orientation along x-axis and G=1"""
    def __init__(self,m1,m2,a):
        self.q = m1/m2
        self.m1 = m1
        self.m2 = m2
        self.M = m1+m2
        self.a = a
        self.x1 = -self.m2/self.M*self.a
        self.x2 = self.m1/self.M*self.a
        self.r1 = np.array([self.x1,0.0,0.0])
        self.r2 = np.array([self.x2,0.0,0.0])
        self.omega = np.sqrt(self.M*self.a**-3)
        self.omega_vec = self.omega * np.array([0,0,1.0])
    
    def phi_roche(self,x,y,z):
        r = np.array([x,y,z])
        phi1 = - self.m1 / np.linalg.norm(r-self.r1)
        phi2 = - self.m2 / np.linalg.norm(r-self.r2)
        cor = - 0.5 * np.linalg.norm( np.cross(self.omega_vec,r)  )**2
        return phi1+phi2+cor
    
    def get_xL_phiL(self,points=1001):
        phi_temp = self.get_phi_function()
        x = np.linspace(-3*self.a,3*self.a,points)
        Lind = argrelextrema(phi_temp(x,0,0),np.greater)[0]
        xL   = x[Lind]
        phiL = phi_temp(x,0,0)[Lind]
        return xL,phiL
    
    def get_phi_function(self):
        return np.vectorize(self.phi_roche)
    
def get_roche_function(orb,time,M1=1,M2=0.3):
    a = np.interp(time,orb['time'],orb['sep'])
    b = makebinary(M1,M2,a)
    xL,phiL = b.get_xL_phiL(points=100)
    return xL, phiL, b.get_phi_function()



def get_plot_array_vertical(quantity,phislice,
                            myfile,profile_file,orb,m1,m2,
                           G=1,rsoft2=0.1,level=0,
                           x1_max=None):
    
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)
    
    
    
    get_cartesian=True
    get_torque=False
    get_energy=False
    if quantity in ['torque_dens_1_z','torque_dens_2_z']:
        get_torque=True
    if quantity in ['ek','ei','etot','epot','epotg','epotp','h','bern']:
        get_energy=True
    
    x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+phislice))]
    d=read_data(myfile,orb,m1,m2,G=G,rsoft2=rsoft2,level=level,
                get_cartesian=get_cartesian,
                get_torque=get_torque,
                get_energy=get_energy,
                x1_max=x1_max,x3_min=x3slicevalue,x3_max=x3slicevalue,
                profile_file=profile_file)
    
    x1 = d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])
    z1 = d['z'][0,:,:]
    val1 = d[quantity][0,:,:]
    
    if(x3slicevalue<0):
        x3slicevalue += np.pi
    else:
        x3slicevalue -= np.pi
    
    d=read_data(myfile,orb,m1,m2,G=G,rsoft2=rsoft2,level=level,
                get_cartesian=get_cartesian,
                get_torque=get_torque,
                get_energy=get_energy,
                x1_max=x1_max,x3_min=x3slicevalue,x3_max=x3slicevalue,
                profile_file=profile_file)
    
    x2 = -d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])
    z2 = d['z'][0,:,:]
    val2 = d[quantity][0,:,:]
    
    # Combine the arrays
    x = np.concatenate((x2,np.flipud(x1)))
    z = np.concatenate((z2,np.flipud(z1)))
    val = np.concatenate((val2,np.flipud(val1)))

    x = np.concatenate((x,x2[0:1]) ,axis=0)
    z = np.concatenate((z,z2[0:1]) ,axis=0)
    val = np.concatenate((val,val2[0:1]) ,axis=0)
    
    return x,z,val

