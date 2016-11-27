# Aleks Diamond-Stanic
# 20161109
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

#setting up arrays with three elements, all zeros - placeholders
wavelengths = [4,8,1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]

xcen = 3629.
ycen = 4153.
dx = 5
dy = 5

# create a PDF file for the plots    
with PdfPages('ad_centroid.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'J0905_final_F4*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(1,3,1)
    plot_image()
    plt.title('F475W')
      
    # read in the F814W image
    file = glob.glob(dir+'J0905_final_F8*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(1,3,2)
    plot_image()
    plt.title('F814W')

    # read in the F160W image
    file = glob.glob(dir+'J0905_final_F1*sci.fits')
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

    os.system('open %s &' % 'ad_centroid.pdf')
    
