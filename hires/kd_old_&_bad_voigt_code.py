#Kwamo 06/20/17
# Attaching Voigt Profile to Galaxy J0905 to Wavelengths of Interest Fe2600 and Mg2796
#Old Voigt profile code, doesn't work 

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.modeling.models import Voigt1D
from scipy.special import wofz
import pylab
from scipy.interpolate import interp1d


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
    #for h in range(0,1):
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
        
        ax = fig.add_subplot(3,1,3)
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        #plt.legend(loc = 1)
        plt.title("Uncertainty of Flux in %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        #handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels)
        
        ax = fig.add_subplot(3,1,2)
        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_tau2796 = column(vel_kms[0], sigma_tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_tau2803 = column(vel_kms[1], sigma_tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')

        g2796 = (vel_2796 > vb) & (vel_2796 < vflip[h])
        g2803 = (vel_2803 > vflip[h]) & (vel_2803 < vr)

        plt.errorbar(vel_2796, col_2796, yerr = sigma_tau2796, linewidth = 0.1, color = '#99ccff', label = 'error')
        plt.errorbar(vel_2803, col_2803, yerr = sigma_tau2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')
        ax.plot(vel_2796, col_2796, color = '#2CA14B', label = names[0])
        ax.plot(vel_2803, col_2803, color = '#2C6EA1', label = names[1])
        plt.title("Uncertainty of Column Density in %s" %(gal[h]))
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Column Density (Particle/cm^2)")
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        plt.legend(loc = 1)


        f = interp1d(vel_kms[0], flux)
        vel_new = np.linspace(-3000, 500, num = 3501, endpoint = True)
        flux_king = f(vel_new)

        voi = Voigt1D(amplitude_L=-0.5, x_0=0.0, fwhm_L=5.0, fwhm_G=5.0)
        #xarr = np.linspace(,5.0,num=40)
        #xarr = vel_kms[0]
        xarr = vel_new
        #yarr = voi(xarr)
        #yarr = flux[good]
        yarr = flux_king - 1.

        voi_init = Voigt1D(amplitude_L=-1.0, x_0=-1000, fwhm_L=200, fwhm_G=200)
        fitter = fitting.LevMarLSQFitter()
        voi_fit = fitter(voi_init, xarr, yarr)
        print(voi_fit)

        ax = fig.add_subplot(3,1,3)
        ax.plot(xarr,yarr+1, color='magenta')
        ax.plot(xarr,voi_fit(xarr)+1, color='red')
        ax.plot(xarr, voi_init(xarr)+1, color='orange')
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        ax.set_ylim (0,2)
        ax.set_xlim(-3000,500)
        
        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("evince  %s &" % 'Error_calculation.pdf')














































##

        
# # read in the relevant SDSS spectrum
# dir = os.environ['SDSSDIR']
# file = dir + 'spec-0761-54524-0409.fits'
# hdulist = fits.open(file)

# # define the coefficients that are used to define the wavelength array
# # For information about how to determine wavelength associated with each pixel: http://classic.sdss.org/dr5/products/spectra/read_spSpec.html
# coeff0 = hdulist[0].header['COEFF0']
# coeff1 = hdulist[0].header['COEFF1']

# # define the flux and wavelength arrays
# flux = hdulist[1].data['flux']
# loglam = hdulist[1].data['loglam']
# model = hdulist[1].data['model']
# ivar = hdulist[1].data['ivar']
# npix = len(flux)
# index = np.arange(npix)

# #To do: check whether this wavelength definition is correct and whether it's in air or vacuum
# #wavelength = 10.**(coeff0 + coeff1*index)

# #Background Information
# wavelength = 10.**loglam
# c = 3 * 10 **8
# z = .7114

# #Background Info Manipulation
# rest_wave = np.array([2600, 2796])
# obs_wave = []
# for i in range(0,len(rest_wave)):
#     wave_calc = (rest_wave[i] * (1+z))
#     obs_wave = np.append(obs_wave,wave_calc)

# velo_int=[]
# for i in range(0,len(obs_wave)):
#     velo_calc = (obs_wave - rest_wave)/(rest_wave) *c
#     velo_int = np.append(velo_int, velo_calc)

# for i in range(0, len(wavelength)):    
#     velocity = np.array([((wavelength*(1+z)) - wavelength[i]) *c / wavelength[i]])





# ## Gaussian function where alpha is the half-width at half-max(Gaussian) and x is the width
# def Gaussian(x, alpha):
#     return np.sqrt(np.log(2)/np.pi)/alpha *np.exp(-(x/alpha)**2*np.log(2))


# ## Lorentzian function where gamma is the half-width at half-max (Lorentzian) and x is the width
# def Lorentzian(x, gamma):
#     return gamma / np.pi / (x**2 + gamma**2)

# #Voigt Function
# def Voigt(x,alpha, gamma):
#     sigma = alpha / np.sqrt(2*np.log(2))
#     return np.real(wofz((x+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi)



# Random Values to see if code works
#alpha, gamma = -200, -200
#x = np.linspace(-0.8,0.8,1000)
#x = np.linspace(-3000,500,3500)
## Gotta find real values for alpha, gamma, and x with Aleks


# graphs Gaussian, Lorentzian, and Voigt
# pylab.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
# pylab.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
# pylab.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
# pylab.xlim(-0.8,0.8)
# pylab.legend()
# pylab.show()

# # create a PDF file for plot output
# filename = 'J0905_sdss.pdf'
# with PdfPages(filename) as pdf:
    
    # fig = plt.figure()


#From ad_voigt Test



##


    
    # ax.plot(xarr,yarr)
    # ax.plot(xarr,voi_fit(xarr))

    # ax = fig.add_subplot(4,1,1)
    # ax.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
    # ax.set_xlim(-3000,500)
    # plt.legend(loc=2)

    # ax = fig.add_subplot(4,1,2)
    # ax.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
    # ax.set_xlim(-3000,500)
    # plt.legend(loc=2)

    # ax = fig.add_subplot(4,1,3)
    # ax.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
    # ax.set_xlim(-3000,500)
    # plt.legend(loc=2)

    # ax = fig.add_subplot(4,1,4)
    # ax.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
    # ax.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
    # ax.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
    # ax.set_xlim(-3000,500)
    # plt.legend(loc=2)

#     fig.tight_layout()
#     # save the figure
#     pdf.savefig()
#     plt.close()
    
#     # open the PDF file
# os.system("open %s &" % filename)

