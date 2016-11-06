# Aleks Diamond-Stanic
# 20161014
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
# dir = os.environ['HSTDIR']
dir = '/Users/sgottlie/data/test/'

wavelengths = [4, 8, 1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]
flux = [0 for x in range(len(wavelengths))] 


# specify the position of the science target
xcen = 3629.
ycen = 4153.

# define the radii to be used for aperture photometry
radii = np.array([1,2,3,4,5,8,10,15,20,25,30,35,40])

# create a PDF file for the plots    
with PdfPages('sg_phot_sed.pdf') as pdf:

    fig = plt.figure()
## HERE BE WHERE I FUCK SHIT UPP

    collection = ['F475W','F814W','F160W']
    # for i in collection:
    for i in range(0, len(collection)):
        
        # read in image
        # how can I modify this code to work for a loop?
        file = glob.glob(dir+'final_' + collection[i] + '*sci.fits')
        print(file)
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']

        # define photometry positions
        positions = [(xcen, ycen)]

        # do photometry
        flux[i] = [0.]
        for radius in radii:
            aperture = CircularAperture(positions, radius)
            phot_table = aperture_photometry(data[i]*fnu[i]/exp[i], aperture)
            flux[i].append(phot_table['aperture_sum'][0])
# HERE BE WHERE I STOP


    # set up the plot
    ax = fig.add_subplot(1,1,1)

    # plot lines connecting photometric points for each aperture
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([flux[0][i+1]-flux[0][i], flux[1][i+1]-flux[1][i], flux[2][i+1]-flux[2][i]])
        ax.plot(wv, fx, color='0.75')

    # plot symbols for each photometric point with color scaled to aperture radius
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([flux[0][i+1]-flux[0][i], flux[1][i+1]-flux[1][i], flux[2][i+1]-flux[2][i]])
        cax = ax.scatter(wv, fx, c=np.array([radii[i],radii[i],radii[i]]), vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2, marker='s')
        
    # set plot parameters
    cbar = fig.colorbar(cax)
    cbar.set_label('radius [pixels]', fontsize=18)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylim([1e-7,1e-4])
    plt.ylabel('flux [Jy]', fontsize=18)
    plt.xlim([300.,2000.])
    plt.xlabel(r'wavelength [$\AA$]', fontsize=18)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'sg_phot_sed.pdf')
