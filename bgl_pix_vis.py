# Aleks Diamond-Stanic
# 20160919
#
# the goals of this code are to do the following:
# (1) read in three images taken in three different filters
# (2) consider the pixels in each image that are in close proximity to the science target
# (3) visualize the colors and flux of these pixels

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
dir = '/Users/adiamond/data/20150203_HST/J0826+4305/coarse/'

# specify the position of the science target and the size of the
# region around the science target to consider
xcen = 3629.
ycen = 4153.
dx = 30.
dy = 30.


# create a PDF file for the plots    
with PdfPages('bgl_pix_vis.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'F475W/final*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp475 = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,1)
    plot_image(stamp475)
    plt.title('F475W')
      
    # read in the F814W image
    file = glob.glob(dir+'F814W/final*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp814 = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,2)
    plot_image(stamp814)
    plt.title('F814W')

    # read in the F160W image
    file = glob.glob(dir+'F160W/final*sci.fits')
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
