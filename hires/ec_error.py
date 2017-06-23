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

# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
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
# define velocity ranges to plot the profile
vb = -3500.
vr = 1500.


minorLocator = AutoMinorLocator()
filename = 'Error calculation.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
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

        #graphs for magnesium absorption lines:
        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_kms[0] > vb) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < vr)
        # plot the profiles using the 2796 profile on the blue side
        # and the 2803 profile on the red side
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0])
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1])
        plt.errorbar(vel_kms[0], flux, yerr = sigma, label = 'error', fmt = 'rs--', ecolor = '#0F0F0F')
        plt.title("Uncertainty")
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 3)
        plt.legend(loc = 1)
        pdf.savefig()
        plt.close()
os.system("open %s &" % 'Error calculation.pdf')
