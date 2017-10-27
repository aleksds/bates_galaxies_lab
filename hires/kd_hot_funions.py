# Kwamae Delva
# 10/27/17
# New Voigt Profile Code Using Linetools and Barak with main goal being to include covering fraction

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy.io import ascii
import sympy
from sympy import mpmath as mp
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting


def voigt_slow(a, u):
      with mp.workdps(20):
            z = mp.mpc(u, a)
            result = mp.exp(-z*z) * mp.erfc(-1j*z)
      return result.real

    # a : float
    #   Ratio of Lorentzian to Gaussian linewidths (see below).
    # u : array of floats, shape (N,)
    #   The frequency or velocity offsets from the line centre, in units
    #   of the FWHM of the Gaussian broadening (see below).
     # workdps(n)(f) :
     #    returns a decorated version of the function f
     #    that sets the decimal precision to n before execution,
     #    and restores the precision afterwards.



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

# define velocity ranges to plot the profile
# vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']
# define velocity ranges to plot the profile


minorLocator = AutoMinorLocator()
filename = 'Tau_Flux_Coldens_Comparison.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']

## Error Bar Calculations

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * np.sqrt(var[i])/fx[i]
            
#Error in Column Density Calculations

        # sigma_coldens = np.zeros([len(flux)])
        # for i in range(0, len(flux)):
        #     sigma_coldens[i] = (sigma[i]/flux[i])

#Velocity Calculations

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5

## Tau Calculations

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        blah = np.log(1/flux)
        tau = blah


## Code to combine 2796 and 2803 absorption profiles

        # 2796 on left (blueshifted side), 2803 on right (redshifted side)
        g2796 = (vel_kms[0] > -3000) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < 500)

        # col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        # col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        # sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        # sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))

        # vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        # vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')

        # for seperated magnesium plots:
        g_2796 = (vel_kms[0] > -3000) & (vel_kms[0] < 500)
        g_2803 = (vel_kms[1] > -3000) & (vel_kms[1] < 500)

        # g2796_coldens = (vel_2796 > -3000) & (vel_2796 < vflip[h])
        # g2803_coldens (vel_2803 > vflip[h]) & (vel_2803 < 500)


        # voigt_slow(vel_kms[0][g2796],flux[g2796])
        p=(2,3,4,5)  #type(q) = tuple
        f=8.123     #type(f) = float
        q=1        #type(q) = integ

                                               ## PLOTTING BEGINS ##
          

                                            ## Seperated Mg II Ion Plots##
                                               
        ## Mg II 2796 


## COmbined Flux Plot

        fig = plt.figure()

        ax = fig.add_subplot(2,1,1)

        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("Flux Plot for Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        # plt.legend(loc = 3)
        
        #Adds error bars to Flux Plot
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)



        fig.tight_layout()
        pdf.savefig()
        plt.close()

        
os.system("open  %s &" % filename)

