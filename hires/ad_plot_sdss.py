# Aleks Diamond-Stanic, 20170606
# code to read in and plot an SDSS spectrum

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
wavelength = 10.**loglam

# create a PDF file for plot output
filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:

    # create a new figure
    fig = plt.figure()

    # make a subplot
    ax = fig.add_subplot(2,1,1)

    # plot flux vs wavelength
    ax.plot(wavelength, flux, linewidth=0.1)
    ax.plot(wavelength, model, linewidth=0.1, color='red')
    ax.plot(wavelength, 1./np.sqrt(ivar), linewidth=0.1)
    
    # set the limits for the x and y axes
    ax.set_ylim(0., 15.)
    ax.set_xlim(3600., 9400)

    # make the axis labels larger
    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.tick_params(axis='y', which='major', labelsize=14)

    # label the x and y axes
    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
    plt.ylabel('$f_{\lambda}$ [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',fontsize=14)

    # save the figure
    pdf.savefig()
    plt.close()
    
    # open the PDF file
    os.system("open %s &" % filename)
