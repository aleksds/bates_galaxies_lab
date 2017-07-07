#Kwamo 06/20/17
# Attaching Voigt Profile to Galaxy J0905 to Wavelengths of Interest Fe2600 and Mg2796

# 20170707 notes from Aleks about things to do:
# (1) figure out how to perform fits with some parameters fixed
# (2) test out having customized parameters for each galaxy (that can be implemented in for loop)
# (3) convert the voigt profile fitting results into column density estimates
# (4) implement fits that are for both Mg II lines simulataneously (with velocity parameters tied together)
# (5) use qualitative and quantitative information from absorption-line profiles and fits to identify trends in the sample

# import relevant packages
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy import constants as const
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.modeling.models import Voigt1D
from scipy.interpolate import interp1d


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
filename = 'voigt_profile_fits.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
    #for h in range(0,1):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * np.sqrt(var[i])/fx[i]

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
        
        ax = fig.add_subplot(3,1,2)
        #plt.errorbar(vel_kms[0], flux, yerr = sigma, label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        #plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)
        ax.plot(vel_kms[0], flux, linewidth=1, label = names[0], color = '#2CA14B')
        #ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        #plt.legend(loc = 1)
        plt.title("Voigt Profile : Velocity vs Flux in %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        #handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels)
        

        f = interp1d(vel_kms[0], flux)
        test = (vel_kms[0] > -3000) & (vel_kms[0] < 500) & (flux < 0.5)
        vel_median = np.median(vel_kms[0][test])
        print(gal[h],vel_median)
        vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)


        blah = len(vel_kms[0][test])
        one = vel_kms[0][test][round(blah*0.2)]
        two = vel_kms[0][test][round(blah*0.4)]
        thr = vel_kms[0][test][round(blah*0.6)]
        fou = vel_kms[0][test][round(blah*0.8)]

        flux_king = f(vel_new)

        xarr = vel_new
        yarr = flux_king - 1.

        voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)
        fitter = fitting.LevMarLSQFitter()
        voi_fit = fitter(voi_init, xarr, yarr)
        print(voi_fit)

        #ax = fig.add_subplot(3,1,3)
        ax.plot(xarr,voi_fit(xarr)+1, color='red')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        ax.set_ylim (0,2)
        ax.set_xlim(-3000,500)
        
        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("evince  %s &" % filename)





























