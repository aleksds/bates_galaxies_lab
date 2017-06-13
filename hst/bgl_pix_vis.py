# Aleks Diamond-Stanic
# 20160919 -- 20161005
#
# Quick description: The code makes "postage stamp" visualizations of
# a galaxy in three different filters and then produces
# color-magntidue diagrams that show the relationship between flux in
# a given filter and the color measured with respect to an adjacent
# filter for all pixels in the postage stamps.
#
# Current status: The current version is set up to work for the galaxy
# J0905.
#
# Future developments: Could tighten up the code with functions and
# for loops.  Could expand to loop over all galaxies.  Could add
# customized apertures.  Could show the location of customized
# apertures on the postage stamp images.
#
# the preliminary goals of this code include the following:
#
# (1) read in three images taken in three different filters
#
# (2) consider the pixels in each image that are in close proximity to
#     the science target
#
# (3) visualize the colors and flux of these pixels
#
#
# more advanced goals include the following:
#
# (A) construct mutltiple apertures (e.g., a central circular aperture
#     and surrounding annuli) and perform photometry in those
#     apertures for all three images
#
# (B) show the photometry from each aperture on each color-magnitude plot
#
# (C) visualize each each aperture on top of the postage stamp images
#
# (D) use the color information in each aperture to estimate a
#     mass-to-light ratio for that aperture


# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
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
xcen = 3388.
ycen = 3504.
dx = 30.
dy = 30.


# create a PDF file for the plots    
with PdfPages('bgl_pix_vis.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'J0905*F4*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp475 = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,1)
    plot_image(stamp475)
    plt.title('F475W')
    plt.suptitle(header['TARGNAME'])
      
    # read in the F814W image
    file = glob.glob(dir+'J0905*F8*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp814 = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,2)
    plot_image(stamp814)
    plt.title('F814W')

    # read in the F160W image
    file = glob.glob(dir+'J0905*F1*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp160 = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,3)
    plot_image(stamp160)
    plt.title('F160W')    

    # plot F160W vs F814W - F160W
    ax = fig.add_subplot(2,3,4)
    ax.plot(-2.5*np.log10(stamp814 / stamp160), np.log10(stamp160+100.)-2., 'ro')
    plt.xlim([-3., 3.])
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.xlabel('-2.5 log (flux[F475W] / flux[F814W])', fontsize=8)
    plt.ylabel('log (flux[F160W]) [counts]', fontsize=8)
    plt.ylim([0.,3.]) 

    # plot F814W vs F475W - F814W
    ax = fig.add_subplot(2,3,6)
    ax.plot(-2.5*np.log10(stamp475 / stamp814), np.log10(stamp814+100.)-2., 'ro')
    plt.xlim([-3., 3.])
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.xlabel('-2.5 log (flux[F475W] / flux[F814W])', fontsize=8)
    plt.ylabel('log (flux[814W]) [counts]', fontsize=8)
    plt.ylim([0.,3.]) 

    ## plot F475W - F814W vs F814W - F160W
    ## goal: clean up the code below, add scaling with flux, add color bar
    #ax = fig.add_subplot(3,3,8)
    #ax.plot(-2.5*np.log10(stamp814 / stamp160), -2.5*np.log10(stamp475 / stamp814), 'ro')
    #plt.xlim([-3., 3.])
    #plt.tick_params(axis='both', which='major', labelsize=8)
    #plt.xlabel('-2.5 log (flux[F814W] / flux[F160W])', fontsize=8)
    #plt.ylabel('-2.5log(F475W/F814W)', fontsize=8)
    #plt.ylim([-3., 3.]) 
    
    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'bgl_pix_vis.pdf')
