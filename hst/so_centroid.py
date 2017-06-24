# Aleks Diamond-Stanic
# 20161109
#
#Senyo Ohene
# Quick description: This code compares centroid values based on three
# routines from photutils: centroid_com, centroid_1dg, centroid_2dg.
#
# Current status: The current version of the code plots these centroid
# values on top of a postage stamp image for the galaxy J0905.
#
# Future developments: Could avoid repetitive code for each filter and
# make a for loop.  Could turn this into a function that returns the
# coordinates of the best centroid location.
#
# import relevant Python modules
import numpy as np
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg

# define a function to plot "postage stamp" images
def plot_image():
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)
    circle0 = plt.Circle((dx,dy),0.1)
    x1, y1 = centroid_com(stamp)
    circle1 = plt.Circle((x1,y1),0.1,color='r')
    ax.add_artist(circle1)
    x2, y2 = centroid_1dg(stamp)
    circle2 = plt.Circle((x2,y2),0.1,color='b')
    ax.add_artist(circle2)
    x3, y3 = centroid_2dg(stamp)
    circle3 = plt.Circle((x3,y3),0.1,color='g')
    ax.add_artist(circle3)
    print(x1, x2, x3)
    print(y1, y2, y3)

# define the directory that contains the images
dir = os.environ['HSTDIR']

# parameters for the initial guess at the galaxy centroid and the size
# of the postage stamps
dx = 5
dy = 5

Galaxies = [['J0826', 3628.7, 4153.8], ['J0901', 3934, 4137], ['J0905', 3386.5, 3503.2], ['J0944', 3477.5, 3404.3], ['J1107', 3572.9, 3339.1], ['J1219', 3803., 4170.], ['J1341', 3885.6, 4164.3], ['J1506', 4149.2, 3921.7], ['J1558', 3787., 4186.], ['J1613', 4175., 3827.], ['J2116', 3566.9, 3435.9], ['J2140', 4067, 4054.4]]
Filters = ['_final_F4*sci.fits', '_final_F8*sci.fits', '_final_F1*sci.fits']
Titles = ['F475W', 'F814W', 'F160W']
# create a PDF file for the plots    
with PdfPages('so_centroid.pdf') as pdf:
    # for each galaxy    
    for j in range(len(Galaxies)):
        fig = plt.figure()
        # for each filter
        for i in range(len(Filters)):
            # read the image
            file = glob.glob(dir+Galaxies[j][0]+Filters[i])
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
    
            # parameters for the initial guess at the galaxy centroid and the size
            xcen = Galaxies[j][1]
            ycen = Galaxies[j][2]
            # create a "postage stamp" image centered on the science target
            stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

            # plot the "postage stamp"
            ax = fig.add_subplot(1,3,i+1)
            plot_image()
            plt.title(Titles[i])
      
         
        plt.suptitle(Galaxies[j][0])
        pdf.savefig()
        plt.close()

os.system('open %s &' % 'so_centroid.pdf')
    

