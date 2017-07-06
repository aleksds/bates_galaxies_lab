#Kwamo 06/20/17
# Attaching Voigt Profile to Galaxy J0905 to Wavelengths of Interest Fe2600 and Mg2796

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.modeling.models import Voigt1D
from scipy.special import wofz
import pylab


# read in the relevant SDSS spectrum
dir = os.environ['SDSSDIR']
file = dir + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)

# define the coefficients that are used to define the wavelength array
# For information about how to determine wavelength associated with each pixel: http://classic.sdss.org/dr5/products/spectra/read_spSpec.html
coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']

# define the flux and wavelength arrays
flux = hdulist[1].data['flux']
loglam = hdulist[1].data['loglam']
model = hdulist[1].data['model']
ivar = hdulist[1].data['ivar']
npix = len(flux)
index = np.arange(npix)

#To do: check whether this wavelength definition is correct and whether it's in air or vacuum
#wavelength = 10.**(coeff0 + coeff1*index)

#Background Information
wavelength = 10.**loglam
c = 3 * 10 **8
z = .7114

#Background Info Manipulation
rest_wave = np.array([2600, 2796])
obs_wave = []
for i in range(0,len(rest_wave)):
    wave_calc = (rest_wave[i] * (1+z))
    obs_wave = np.append(obs_wave,wave_calc)

velo_int=[]
for i in range(0,len(obs_wave)):
    velo_calc = (obs_wave - rest_wave)/(rest_wave) *c
    velo_int = np.append(velo_int, velo_calc)

for i in range(0, len(wavelength)):    
    velocity = np.array([((wavelength*(1+z)) - wavelength[i]) *c / wavelength[i]])



## Gaussian function where alpha is the half-width at half-max(Gaussian) and x is the width
def Gaussian(x, alpha):
    return np.sqrt(np.log(2)/np.pi)/alpha *np.exp(-(x/alpha)**2*np.log(2))


## Lorentzian function where gamma is the half-width at half-max (Lorentzian) and x is the width
def Lorentzian(x, gamma):
    return gamma / np.pi / (x**2 + gamma**2)

#Voigt Function
def Voigt(x,alpha, gamma):
    sigma = alpha / np.sqrt(2*np.log(2))
    return np.real(wofz((x+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi)



# Random Values to see if code works
alpha, gamma = 200, 200
#x = np.linspace(-0.8,0.8,1000)
x = np.linspace(-3000,500,3500)
## Gotta find real values for alpha, gamma, and x with Aleks


# graphs Gaussian, Lorentzian, and Voigt
# pylab.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
# pylab.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
# pylab.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
# pylab.xlim(-0.8,0.8)
# pylab.legend()
# pylab.show()

# create a PDF file for plot output
filename = 'J0905_sdss.pdf'
with PdfPages(filename) as pdf:
    
    fig = plt.figure()
    
    ax = fig.add_subplot(4,1,1)
    ax.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
    ax.set_xlim(-3000,500)
    plt.legend(loc=2)

    ax = fig.add_subplot(4,1,2)
    ax.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
    ax.set_xlim(-3000,500)
    plt.legend(loc=2)

    ax = fig.add_subplot(4,1,3)
    ax.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
    ax.set_xlim(-3000,500)
    plt.legend(loc=2)

    ax = fig.add_subplot(4,1,4)
    ax.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
    ax.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
    ax.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
    ax.set_xlim(-3000,500)
    plt.legend(loc=2)

    fig.tight_layout()
    # save the figure
    pdf.savefig()
    plt.close()
    
    # open the PDF file
os.system("open %s &" % filename)

