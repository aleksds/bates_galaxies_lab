# Aleks Diamond-Stanic
# 20161109
#
# Quick description: This code compares centroid values based on three
# routines from photutils: centroid_com, centroid_1dg, centroid_2dg.
#
# Current status: The current version of the code plots these centroid
# values on top of a postage stamp image for the galaxy J0905.
#
# Future developments: Could turn this into a function that returns the
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
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    return [((x1-dx+xcen),(x2-dx+xcen),(x3-dx+xcen)) , ((y1-dy+ycen),(y2-dy+ycen),(y3-dy+ycen))]

# define the directory that contains the images
dir = os.environ['HSTDIR']

# parameters for the initial guess at the galaxy centroid and the size
Galaxies = [['J0826', 3628.7, 4153.8], ['J0901', 3934, 4137], ['J0905', 3386.5, 3503.2], ['J0944', 3477.5, 3404.3], ['J1107', 3572.9, 3339.1], ['J1219', 3803., 4170.], ['J1341', 3885.6, 4164.3], ['J1506', 4149.2, 3921.7], ['J1558', 3787., 4186.], ['J1613', 4175., 3827.], ['J2116', 3566.9, 3435.9], ['J2140', 4067, 4054.4]]
xcen = 3386.5
ycen = 3503.2
dx = 5
dy = 5

# create a PDF file for the plots
foureightone = [4,8,1]
collection = ['F475W','F814W','F160W']
with PdfPages('ad_centroid.pdf') as pdf:

    fig = plt.figure()
    plt.scatter(xcen,ycen, label='current')
    yolo = ['centroid_com','centroid_1dg','centroid_2dg']
    for i in range(0,len(collection)):
        file = glob.glob(dir+'J0905_final_F'+str(foureightone[i])+'*sci.fits')
        hdu = fits.open(file[0])
        data, header = hdu[0].data, hdu[0].header
        stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]
        coor = plot_image()
        plt.scatter(coor[0],coor[1],label = yolo[i])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(3385,3387)
    plt.ylim(3502.5,3504)
    plt.title('Centroid Calculations')
    legend = plt.legend(loc='upper left', shadow=True)
        
    pdf.savefig()
    plt.close()

    os.system('open %s &' % 'ad_centroid.pdf')
