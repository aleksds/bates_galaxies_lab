# Aleks Diamond-Stanic
# 20160916
#
# Quick description: This code reads in images of the same galaxy in
# three filters and makes "postage stamp" visualizations with the goal
# of providing a qualitative sense for how the morphology of the
# galaxy and its luminosity change as a function of wavelength. This
# code is expanded upon in bgl_aper_phot.py and bgl_pix_vis.py.
#
# Current status: The current version is set up to work for the galaxy
# J0826.
#
# Future developments: Could tighten up the code with functions and
# for loops.  Could expand to loop over all galaxies.  Could produce
# three-color images.  
#
# Original notes:
#
# the goals of this code are to do the following:
# (1) read in three images taken in three different filters
# (2) plot a "postage stamp" version of each image, centered on the
#     science target
# (3) combine the three image to make a single RGB image (not yet implemented)
#
#     For future reference, here are some relevant links (from a
#     Google search for  "astropy make a three color image"):
#
#     (A) "Making RGB images from FITS files with python / matplotlib"
#     http://www.astrobetter.com/blog/2010/10/22/making-rgb-images-from-fits-files-with-pythonmatplotlib/
#
#     (B) "Making a publication quality image"
#     https://python4astronomers.github.io/intro/quick-tour.html
#
# import relevant Python modules
import numpy as np
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# define a function to plot "postage stamp" images
def plot_image():
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)

# define the directory that contains the images
dir = os.environ['HSTDIR']
xcen = 3629.
ycen = 4153.
dx = 500
dy = 500

# create a PDF file for the plots    
with PdfPages('bgl_image_stamp.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'J0826*F475W*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(1,3,1)
    plot_image()
    plt.title('F475W')
      
    # read in the F814W image
    file = glob.glob(dir+'J0826*F814W*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(1,3,2)
    plot_image()
    plt.title('F814W')

    # read in the F160W image
    file = glob.glob(dir+'J0826*F160W*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(1,3,3)
    plot_image()
    plt.title('F160W')    
    
    pdf.savefig()
    plt.close()

os.system('open %s &' % 'bgl_image_stamp.pdf')
    
