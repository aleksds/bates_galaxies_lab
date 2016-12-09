# Aleks Diamond-Stanic
# 20160920 -- 20161005
#
# the preliminary goals of this code are to do the following:
#
# (1) read in three images taken in three different filters
#
# (2) consider the pixels in each image that are in close proximity to
#     the science target
#
# (3) perform photometry in circular apertures
#
#
# the more advanced goals include the following:
#
# (A) query the SDSS database for ugriz photometry for the science
#     target
#
# (B) use the SDSS photometry to begin to construct an SED that will
#     extend between 0.2 and 2 microns in the observed frame
#
# (C) show the profile of each filter bandpass on the SED plot
#
# (D) measure the total photometry in the F475W, F814W, and F160W
#     filters for the scence target
#
# (E) include the three HST photometry points on the broad-band SED
#
# (F) construct a nuclear aperture, an intermediate annulus, and an
#     outer annulus for the HST images
#
# (G) measure photometry in each of these three apertures for each of
#     the three HST filters
#
# (H) construct three-band SEDs for each aperture and include these
#     SEDs on the plot

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

# define a function to plot "postage stamp" images
def plot_image(stamp):
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)

# define the directory that contains the images
dir = os.environ['HSTDIR']

# specify the position of the science target and the size of the
# region around the science target to consider
xcen = 3629.
ycen = 4153.
dx = 100.
dy = 100.


# create a PDF file for the plots    
with PdfPages('bgl_aper_phot.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'J0826*F475W*sci.fits')
    hdu = fits.open(file[0])
    data475, header475 = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp475 = data475[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,1)
    plot_image(stamp475)
    plt.title('F475W')
    plt.suptitle(header475['TARGNAME'])
      
    # read in the F814W imagex
    file = glob.glob(dir+'J0826*F814W*sci.fits')
    hdu = fits.open(file[0])
    data814, header814 = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp814 = data814[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,2)
    plot_image(stamp814)
    plt.title('F814W')

    # read in the F160W image
    file = glob.glob(dir+'J0826*F160W*sci.fits')
    hdu = fits.open(file[0])
    data160, header160 = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp160 = data160[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,3)
    plot_image(stamp160)
    plt.title('F160W')    

    # plot F160W photometric curve of growth
    positions = [(xcen, ycen)]
    radii = np.arange(dx)+1
    flux = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data160, aperture)
        flux.append(phot_table['aperture_sum'][0])
    ax = fig.add_subplot(2,3,6)
    ax.plot(radii, flux, 'ro')
    plt.ylim([0,4e6])
    plt.xlabel('aperture radius')

    # plot F814W photometric curve of growth
    positions = [(xcen, ycen)]
    radii = np.arange(dx)+1
    flux = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data814, aperture)
        flux.append(phot_table['aperture_sum'][0])
    ax = fig.add_subplot(2,3,5)
    ax.plot(radii, flux, 'ro')
    plt.ylim([0,4e6])
    plt.xlabel('aperture radius')
    
    # plot F475W photometric curve of growth
    positions = [(xcen, ycen)]
    radii = np.arange(dx)+1
    flux = []
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data475, aperture)
        flux.append(phot_table['aperture_sum'][0])
    ax = fig.add_subplot(2,3,4)
    ax.plot(radii, flux, 'ro')
    plt.ylim([0,4e6])
    plt.xlabel('aperture radius')
    
    pdf.savefig()
    plt.close()

    
os.system('open %s &' % 'bgl_aper_phot.pdf')
