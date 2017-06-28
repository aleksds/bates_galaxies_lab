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

def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
    return cor_col

# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = 3E5
#mass of electron
mass_e = 9.11**(-28)

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
# define velocity ranges to plot the profile
vb = -3000.
vr = 500.


minorLocator = AutoMinorLocator()
filename = 'Error_calculation.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
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
        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        blah = np.log(1/flux)
        tau = blah
        sigma_tau = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma_tau[i] = (sigma[i]/flux[i])
                
        #graphs for magnesium absorption lines:
        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_kms[0] > vb) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < vr)
        # plot the profiles using the 2796 profile on the blue side
        # and the 2803 profile on the red side
        fig = plt.figure()
        
        ax = fig.add_subplot(2,1,1)
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#8BFEFA', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#8BFEFA', markevery = 10, linewidth = .1)
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        #plt.legend(loc = 1)
        plt.title("Uncertainty of Flux in %s" %(gal[h]))
        plt.xlabel("Velocity")
        plt.ylabel("Flux")
        #handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels)
        
        ax = fig.add_subplot(2,1,2)
        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_tau2796 = column(vel_kms[0], sigma_tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_tau2803 = column(vel_kms[1], sigma_tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')

        g2796 = (vel_2796 > vb) & (vel_2796 < vflip[h])
        g2803 = (vel_2803 > vflip[h]) & (vel_2803 < vr)

        plt.errorbar(vel_2796, col_2796, yerr = sigma_tau2796, linewidth = 0.1, color = '#8BFEFA', label = 'error')
        plt.errorbar(vel_2803, col_2803, yerr = sigma_tau2803, linewidth = 0.1, color = '#8BFEFA', label = '_nolegend_')
        ax.plot(vel_2796, col_2796, color = '#2CA14B', label = names[0])
        ax.plot(vel_2803, col_2803, color = '#2C6EA1', label = names[1])
        plt.title("Uncertainty of Column Density in %s" %(gal[h]))
        plt.xlabel("Velocity")
        plt.ylabel("Column Density")
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        plt.legend(loc = 1)
        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("evince  %s &" % 'Error_calculation.pdf')
