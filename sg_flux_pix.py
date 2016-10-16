# Sophia CW Gottlieb I
# 20161014
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
dir = '/Users/sgottlie/data/test/'

# specify the position of the science target and the size of the
# region around the science target to consider
xcen = 3629.
ycen = 4153.
dx = 100.
dy = 100.

# create a PDF file for the plots    
with PdfPages('sg_flux_pix.pdf') as pdf:
    fig = plt.figure()

    
    collection = ['F475W','F814W','F160W']
    filter = [475, 814, 1600]
    colors = ['bo', 'go', 'ro']
    
    for i in range(0, len(collection)):
        
        # read in image
        # how can I modify this code to work for a loop?
        file = glob.glob(dir+'final_' + collection[i] + '*sci.fits')
        print(file)
        hdu = fits.open(file[0])
        data, header = hdu[0].data, hdu[0].header

        plt.suptitle(header['TARGNAME'])

        #plot the thing; do the stuff
        positions = [(xcen, ycen)]
        radii = np.arange(dx)+1
        wavelength = np.arange(dx)+1
        areaflux = np.arange(dx)+1

        flux = []
        for radius in radii:
            aperture = CircularAperture(positions, radius)
            phot_table = aperture_photometry(data, aperture)
            flux.append(phot_table['aperture_sum'][0])
        cx = fig.add_subplot(3,1,1)
        ax = fig.add_subplot(3,1,2)
        bx = fig.add_subplot(3,1,3)

        wavelength[0] = filter[i]
        areaflux[0] = flux[0]
        for j in range (1, len(radii)):
            wavelength[j] = filter[i]
            areaflux[j] = flux[j] - flux [j-1]
        cx.plot(wavelength, np.log10(areaflux), colors[i])
        # ax.plot(wavelength, flux, 'ro')
        # plt.ylim([0,4e6])
        #plt.ylim([0,6])
        plt.xlim([400,1675])
        plt.xlabel('wavelength in nm')
        plt.ylabel('mother fluxer')
        plt.tick_params(axis='both',which='major', labelsize = 8)
        fig.tight_layout()
                
        #do stuff
        ax.plot(radii, flux, colors[i])
        plt.xlabel('aperture radius')
        plt.ylabel('cumulative flux')
        plt.tick_params(axis='both',which='major', labelsize = 8)
        fig.tight_layout()
        
        bx.plot(radii, areaflux, colors[i])
        plt.xlabel('aperture radius')
        plt.ylabel('subtractive flux')
        plt.xlim([0,100])
        plt.tick_params(axis='both',which='major', labelsize = 8)
        fig.tight_layout()
        

    pdf.savefig()
    plt.close()
os.system('open %s &' % 'sg_flux_pix.pdf')

