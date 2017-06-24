#Kwamo 06/20/17
# Attaching Voigt Profile to Galaxy J0905 to Wavelengths of Interest Fe2600 and Mg2796
Chill
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
    return np.sqrt(np.log(2)/np.pi)/alpha *np.exp(-(x/alpha**2*np.log(2)))


## Lorentzian function where gamma is the half-width at half-max (Lorentzian) and x is the width
def Lorentzian(x, gamma):
    return gamma / np.pi / (x**2 + gamma**2)

#Voigt Function
def Voigt(x,alpha, gamma):
    sigma = alpha / np.sqrt(2*np.log(2))
    return np.real(wofz((x+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi)



# alpha for Gaussian Function    Formula is sigma * sqrt (2* ln(2) ) ,  where sigma is std dev
alpha = np.std(flux) * np.sqrt(2*np.log(2))

# gamma for Lorentzian Function    Formula is 
#gamma = 

# x accounts for the width of the curve based on error/deviations
# for i in velo_int:
#     x = ((np.std(velocity) *3 + velo_int)  - velo_int) #(Ask Aleks about frequency in this case)

x = velocity



# # create a PDF file for plot output
# filename = 'J0905_sdss.pdf'
# with PdfPages(filename) as pdf:

# #     # fig = plt.figure
# #     pylab.plot(x, Gaussian(x, alpha), ls=':', c='k', label='Gaussian')
# #     pylab.plot(x, Lorentzian(x, gamma), ls='--', c='k', label='Lorentzian')
# #     pylab.plot(x, Voigt(x, alpha, gamma), c='k', label='Voigt')
# #     pylab.xlim(np.min(velocity),np.max(velocity))
# #     pylab.legend()
# #     pylab.show()
#     # ## Gaussian Wavelength

#     # fig = plt.figure()
#     # ax = fig.add_subplot (1,1,1)
#     # ax.plot(wavelength,g_wave(wavelength), label = 'Gaussian Fit of Wavelength')
#     # plt.xlabel('Wavelength')
#     # plt.ylabel('Flux')
#     # plt.legend(loc=2)

#     # # save the figure
#     # pdf.savefig()
#     # plt.close()

#     # ##Gaussian Velocity
    
#     # fig = plt.figure()
#     # ax = fig.add_subplot (1,1,1)
#     # ax.plot(velocity,g_velo(velocity), label = 'Gaussian Fit of Velocity')
#     # plt.xlabel('Velocity')
#     # plt.ylabel('Flux')
#     # plt.legend(loc=2)

#     # save the figure
#     pdf.savefig()
#     plt.close()
    
#     # open the PDF file
#     os.system("open %s &" % filename)

