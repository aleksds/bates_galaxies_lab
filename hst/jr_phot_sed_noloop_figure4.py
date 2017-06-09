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
dir = os.environ['HSTDIR']

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3388.
ycen = 3504.
dx = 100
dy = 100

# define the radii to be used for aperture photometry
radii = np.array([1,10,20,30,40])

# create a PDF file for the plots    
with PdfPages('jr_phot_sed_noloop.pdf') as pdf:
    
    fig = plt.figure()
    
    # read in the F475W image
    file = glob.glob(dir+'final_F4*sci.fits')
    hdu = fits.open(file[0])
    data475, header475 = hdu[0].data, hdu[0].header

    # read in the F814W image
    file = glob.glob(dir+'final_F8*sci.fits')
    hdu = fits.open(file[0])
    data814, header814 = hdu[0].data, hdu[0].header

    # read in the F160W image
    file = glob.glob(dir+'final_F1*sci.fits')
    hdu = fits.open(file[0])
    data160, header160 = hdu[0].data, hdu[0].header

    # define positions for photometry
    positions = [(xcen, ycen)]
    
    # do photometry on F475W image
    flux475 = [0.]
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data475, aperture)
        flux475.append(phot_table['aperture_sum'][0])
        
    # do photometry on F814W image
    flux814 = [0.]
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data814, aperture)
        flux814.append(phot_table['aperture_sum'][0])
        
    # do photometry on F160W image
    flux160 = [0.]
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data160, aperture)
        flux160.append(phot_table['aperture_sum'][0])

    # set up the plot
    ax = fig.add_subplot(1,1,1)

    # plot lines connecting photometric points for each aperture
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([flux475[i+1]-flux475[i], flux814[i+1]-flux814[i], flux160[i+1]-flux160[i]])
        ax.plot(wv, fx, color='0.75')

    # plot symbols for each photometric point with color scaled to aperture radius
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([flux475[i+1]-flux475[i], flux814[i+1]-flux814[i], flux160[i+1]-flux160[i]])
        cax = ax.scatter(wv, fx, c=np.array([radii[i],radii[i],radii[i]]), vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2, marker='s')
        plt.title('Flux vs. Wavelength, Broad Annuli')
        
    # set plot parameters
    cbar = fig.colorbar(cax)
    cbar.set_label('radius [pixels]', fontsize=18)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylim([1e3,6e7])
    plt.ylabel('flux [image units]', fontsize=18)
    plt.xlim([300.,2000.])
    plt.xlabel(r'wavelength [$\AA$]', fontsize=18)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_phot_sed_noloop.pdf')
