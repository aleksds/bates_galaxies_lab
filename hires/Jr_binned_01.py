#Jose
#Binned plot for flux and velocity

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator


# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('km/s')

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA

#lines of interest
lines=[mgii2796,mgii2803]
# arrays of galaxy names, redshifts, approximate centroid velocities,
# approximate maximum velocities, and velocities to "flip" between the 2796
# and 2803 line when constructing an empirical mgii velocity profile
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip'] * u.km / u.s
lines=[mgii2803,mgii2796]
# sort by maximum outflow velocity
vmax_index = np.flipud(np.argsort(vmax))




# setting parameters
minv=-3500#* u.km / u.s #km/s
maxv=500#* u.km / u.s #km/s
deltav=30#* u.km / u.s #km/s binning

#Make plot for binned velocity and flux

filename='Binned.vel_flux'+str(deltav)+'.pdf'
xls=5
yls=5
minorLocator=AutoMinorLocator()
with PdfPages(filename) as pdf:
    for i in range(0,len(gal)):

        #read in the spectrum
        datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        #Intensity=1-flux



        #create velocity array
        grid_vel=np.arange(minv,maxv,deltav)
        #create storing arrays
        med_vel=np.zeros(len(grid_vel))
        flux_vel=np.zeros(len(grid_vel))
        std_vel=np.zeros(len(grid_vel))
        vel=np.zeros([len(lines),len(wave)])

        for k in range(0,len(lines)):
            vel[k]=((wave-lines[k]*(1+zem[i]))/(lines[k]*(1+zem[i])))*c
      

#maybe I have to run different loops for med,mean,and std...? 
            for j in range(0,len(grid_vel)):
            
                #Identify elements in desired velocity range
                good=((vel[k]>grid_vel[j]-10) & (vel[k]<grid_vel[j]+10))
                med_vel[j]=np.median(flux[good])
                flux_vel[j]=np.mean(flux[good])
                std_vel[j]=np.std(flux[good])
                #print(j, grid_vel[j], np.sum(good))


                
        fig = plt.figure()
                
        plt.scatter(grid_vel,flux_vel) #still need to figure out how to plot error bars
        plt.errorbar(grid_vel,flux_vel,yerr=std_vel)
        plt.xlim([minv,maxv])
        
        pdf.savefig()
        plt.close()

os.system("open %s &" % filename)
