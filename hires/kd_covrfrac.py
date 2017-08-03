# Makes plot for covering fraction vs velocity

from scipy.interpolate import spline
import math
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch


# define the data directory
dir = os.environ['HIRESDIR']

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853
mgii2796 = 2796.3542699
#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
# array of names, lines of interest, oscillator strength:
names = ['Mg II 2796', 'Mg II 2803']
lines = [mgii2796, mgii2803]
fosc = [f2796, f2803]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']


minorLocator = AutoMinorLocator()
filename = 'Cover_Fraction.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
    # for h in range(0, 1):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5
    
        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_kms[0] > -3000) & (vel_kms[0] <500)#< vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < 500)
        # plots the profiles where the 2796 profile on the blue side
        # and the 2803 profile on the red side

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        #Used 2803 line because all graphs have their max absorption in Mg2803 wavelength region
        f_new_2796 = ((flux - np.min(flux[g2796]))/(1 - np.min(flux[g2796])))
        f_new_2803 = ((flux - np.min(flux[g2803]))/(1 - np.min(flux[g2803])))
        # blah = np.log(1/flux)
        bleh = np.log (1/f_new_2796)
        blah = np.log (1/f_new_2803)
        tau_2796 = bleh
        tau_2803 = blah

        fig = plt.figure()

##Flux Plots
        
        ax = fig.add_subplot(2,1,1)
        plt.title("Flux Plot in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)


## Tau Plots

        ax = fig.add_subplot(2,1,2)
        plt.title("Tau based on constant Covering Fraction in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Tau")
        ax.plot(vel_kms[0][g2796], tau_2796[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], tau_2803[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        # y-axis limit set to 5 because tau of 5 and t of infinity have no visible difference
        ax.set_ylim(0, 2)


        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("open %s &" % 'Cover_Fraction.pdf')
