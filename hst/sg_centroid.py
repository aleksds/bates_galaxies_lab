# Sophia C. W. Gottlieb I
# 20170208
#
# Adapted from
# Aleks Diamond-Stanic
# 20161109
#
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
# PLOT IMAGE GIVES VALUES OF PIXELS? 4-5? WHEN COMPOVERLAY NEEDS 3000-4000??
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
    return ((x1, y1),(x2, y2),(x3, y3))
    #print(x1, x2, x3)
    #print('now some y values')
    #print(y1, y2, y3)

# define the directory that contains the images
dir = os.environ['HSTDIR']

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
w = 2

#I THINK THESE NEED TO BE CHANGED
# parameters for the initial guess at the galaxy centroid and the size
# of the postage stamps
#add 
xcens = [3629, 3933, 3386, 3477, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycens = [4154, 4136, 3503, 3404, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]

xadds = [0,1, 0, 0, 0,2.0,-2.2,-1.5,-0.5, 1.7,1,0]
yadds = [0,0,-1,-2, 1,1.5, 1.0, 1.5,-0.6, 1  ,1,0]
dx = 5
dy = 5
collection = ['F475W','F814W','F160W']
with PdfPages('sg_centroid.pdf') as pdf:
    # create a PDF file for the plots    
    for w in range(0,len(galaxies)):
    #for w in range(11,12):
        fig = plt.figure()
        plt.suptitle(galaxies[w])
        xcen = xcens[w]+xadds[w]
        ycen = ycens[w]+yadds[w]
        
        for i in range(0,len(collection)):
            # read in the F814W image
            file = glob.glob(dir+galaxies[w]+'_final_'+collection[i]+'*sci.fits')
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
    
            # create a "postage stamp" image centered on the science target
            stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

            # plot the "postage stamp"
            ax = fig.add_subplot(1,3,i+1)
            plot_image()
            plt.title(collection[i])    
    
        pdf.savefig()
        plt.close()

    os.system('open %s &' % 'sg_centroid.pdf')
    
