#Kwamae Delva 6/20/17
#Shows the absorption profiles for Mg 2796 and Mg 2803

from scipy.interpolate import spline
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator


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
filename = 'MgII_Absorption_Profile.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
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
        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        blah = np.log(1/flux)
        tau = blah
                
        #graphs for magnesium absorption lines:
        g_2796 = (vel_kms[0] > -3000) & (vel_kms[0] < 500)
        g_2803 = (vel_kms[1] > -3000) & (vel_kms[1] < 500)

        fig = plt.figure()
        
        ax = fig.add_subplot(2,1,1)
        ax.plot(vel_kms[0][g_2796], flux[g_2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("MgII 2796 Flux Plot in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")


        ax = fig.add_subplot(2,1,2)
        ax.plot(vel_kms[1][g_2803], flux[g_2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("MgII 2803 Flux Plot in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")

        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("open  %s &" % 'MgII_Absorption_Profile.pdf')
