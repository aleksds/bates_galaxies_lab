# Aleks Diamond-Stanic
# 20161014
#
# the main goal of this code is to perform photometry in different
# apertures in three different images for a galaxy and then plots flux
# vs wavelegnth (i.e., a spectral energy distribution) for each aperture

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# define the directory that contains the images
dir = '/Users/adiamond/data/20150203_HST/J0826+4305/coarse/'

# specify the position of the science target and the size of the
# region around the science target to consider
xcen = 3629.
ycen = 4153.
#dx = 100.
#dy = 100.
radii = np.array([1,2,3,4,5,8,10,15,20,25,30,35,40])

# create a PDF file for the plots    
with PdfPages('ad_phot_sed.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'F475W/final*sci.fits')
    hdu = fits.open(file[0])
    data475, header475 = hdu[0].data, hdu[0].header

    # read in the F814W imagex
    file = glob.glob(dir+'F814W/final*sci.fits')
    hdu = fits.open(file[0])
    data814, header814 = hdu[0].data, hdu[0].header

    # read in the F160W image
    file = glob.glob(dir+'F160W/final*sci.fits')
    hdu = fits.open(file[0])
    data160, header160 = hdu[0].data, hdu[0].header

    ax = fig.add_subplot(1,1,1)
    positions = [(xcen, ycen)]
    # plot F475W photometric points
    flux475 = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data475, aperture)
        flux475.append(phot_table['aperture_sum'][0])
        wv = np.array([0.,475])
        fx = np.array([0., flux475[-1]])
        ax.plot(wv, fx, 'ro')

    # plot F814W photometric points
    flux814 = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data814, aperture)
        flux814.append(phot_table['aperture_sum'][0])
        wv = np.array([0.,814])
        fx = np.array([0., flux814[-1]])
        ax.plot(wv, fx, 'ro')
        
    # plot F160W photometric points
    flux160 = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data160, aperture)
        flux160.append(phot_table['aperture_sum'][0])
        wv = np.array([0.,1600])
        fx = np.array([0., flux160[-1]])
        ax.plot(wv, fx, 'ro')    
        
    # set plot parameters
    plt.ylim([0,4e6])
    plt.xlim([300.,2000.])
    plt.xlabel('wavelength')  

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'ad_phot_sed.pdf')
