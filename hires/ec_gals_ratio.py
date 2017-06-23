#This program loops through all the galaxies we are focused on
#2796 and 2803 Mg 
#2600 and 2382 Fe
#Creates a page with four graphs for each galaxy

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
import random
import math


# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

# define velocity ranges to plot the profile
vb = -3500.
vr = 1500.

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
feii2600 = 2600.1724835 * u.AA
feii2382 = 2382.7641781 * u.AA
 
#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
f2600 = 0.2394
f2382 = 0.320

# array of names, lines of interest, oscillator strength:
names = ['Mg II 2796', 'Mg II 2803', 'Fe II 2600','Fe II 2382']
lines = [mgii2796, mgii2803, feii2600, feii2382]
fosc = [f2796, f2803, f2600, f2382]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip'] 

def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
    return cor_col

def velocity(vel, col_dens):
    cor_vel = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_vel = np.append(cor_vel, vel[i])
    return cor_vel

counter1 = 0
minorLocator = AutoMinorLocator()
filename = 'Ratio.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
        counter1 += 1
        counter = 0
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        fx = data['fx']
        var = data['var']

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * math.sqrt(var[i])/fx[i]

        
        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * c_kms
        tau = np.zeros([len(lines),len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        for j in range(0, len(lines)):
            blah = np.log(1/flux)
            tau[j] = blah
        #graphs for the iron absorption lines(flux):    
        col_dens = np.zeros([len(lines), len(flux)])
        #calculating column density
        fig = plt.figure()
        
        ax = fig.add_subplot(1,1,1)
        for k in range(2, 4):
            f = interp1d(vel_kms[k], flux)
            vel_new = np.linspace(-3000, 500, num = 3501, endpoint = True)
            flux_king = f(vel_new)
            ax.plot(vel_new, flux_king, lw = .9, label = names[k])
            
            plt.ylabel('Flux')
            ax.set_xlim(-3000, 500)
            ax.set_ylim(0, 2)
            plt.title("Flux of %s" % gal[h])

        #graphs for magnesium absorption lines:
        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_kms[0] > vb) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < vr)
        # plot the profiles using the 2796 profile on the blue side
        # and the 2803 profile on the red side
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0])
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1])
        #plt.text(xmin+0.03*(xmax-xmin), 0.15, gal[indx])
        plt.legend(loc = 1)
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.errorbar(vel_kms[0], flux, yerr = sigma, label = 'error', fmt = 'rs--')
        plt.title("Uncertainty")
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 1)
        pdf.savefig()
        plt.close()
        
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #column density graphs
        for k in range(0, 4):
            col_dens[k] = tau[k] / (2.654E-15 * fosc[k] * (wave/(1+zem[h])) * fosc[k])
            cor_col = column(vel_kms[k], col_dens[k])
            #f = interp1d(vel_kms[k], col_dens[k])
            vel_new = np.linspace(-3000, 500, num = len(cor_col), endpoint = True)
            #col_d = f(vel_new)
            ax.plot(vel_new, cor_col, lw = .9, label = names[k])
            plt.ylabel('Column Density')
            ax.set_xlim(-3000, 500)
            #ax.set_ylim(0, 2)
            plt.title("Column Density of %s" % gal[h])
            
        plt.legend(loc = 1)
        pdf.savefig()
        plt.close()
            
os.system("open %s &" % filename)



