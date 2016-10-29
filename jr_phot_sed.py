# Joshua Rines
# 20161029
#
# the main goal of this code is to perform photometry in different
# apertures in three different images for a galaxy and then plots flux
# vs wavelegnth (i.e., a spectral energy distribution) for each aperture

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

# define the directory that contains the images
dir = '/Users/jrines/data/test/'

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3388.
ycen = 3504.
dx = 100
dy = 100

# create a PDF file for the plots    
with PdfPages('jr_phot_sed.pdf') as pdf:

    fig = plt.figure()

    collection = ['F470W','F814W','F160W']
    filter = [475,814,1600]

    for i in range (0, len(collection)):

        # define the radii to be used for aperture photometry
        radii = np.array([1,2,3,4,5,8,10,15,20,25,30,35,40])
        
        # read in the images
        file = glob.glob(dir+'final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data475, header475 = hdu[0].data, hdu[0].header

        #define positions for photometry
        positions = [(xcen, ycen)]

        #do photometry on images
        for j in range (0, len(filter)):
            flux[j] = [0.]
            for radius in radii:
                aperature = CircularAperture(positions, radius)
                phot_table = aperture_photometry(data[j], aperture)
                flux[j].append(phot_table['aperture_sum'][0])
        
        #set up the plot
        ax = fig.add_subplot(1,1,1)

        # plot lines connecting photometric points for each aperture
        for k in range(0,len(radii)):
            wv = np.array([475.,814.,1600])
            fx = np.array([flux475[k+1]-flux475[k], flux814[k+1]-flux814[k], flux160[k+1]-flux160[k]])
            ax.plot(wv, fx, color='0.75')

        # plot symbols for each photometric point with color scaled to aperture radius
        for i in range(0,len(radii)):
            wv = np.array([475.,814.,1600])
            fx = np.array([flux475[i+1]-flux475[i], flux814[i+1]-flux814[i], flux160[i+1]-flux160[i]])
            cax = ax.scatter(wv, fx, c=np.array([radii[i],radii[i],radii[i]]), vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2, marker='s')

        # set plot parameters
        cbar = fig.colorbar(cax)
        cbar.set_label('radius [pixels]', fontsize=18)
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.ylim([1e4,6e5])
        plt.ylabel('flux [image units]', fontsize=18)
        plt.xlim([300.,2000.])
        plt.xlabel(r'wavelength [$\AA$]', fontsize=18)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_phot_sed.pdf')
