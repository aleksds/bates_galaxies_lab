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
dir = os.environ['HSTDIR']

# specify the position of the science target
xcen = 3629.
ycen = 4153.

# define the radii to be used for aperture photometry
radii = np.array([1,2,3,4,5,8,10,15,20,25,30,35,40])

# create a PDF file for the plots    
with PdfPages('ad_phot_sed.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'J0826*F475W*sci.fits')
    hdu = fits.open(file[0])
    data475, header475 = hdu[0].data, hdu[0].header

    # read in the F814W image
    file = glob.glob(dir+'J0826*F814W*sci.fits')
    hdu = fits.open(file[0])
    data814, header814 = hdu[0].data, hdu[0].header

    # read in the F160W image
    file = glob.glob(dir+'J0826*F160W*sci.fits')
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

    photfnu475 = header475['PHOTFNU']
    exptime475 = header475['EXPTIME']
    fnujy475 = np.array(flux475) * photfnu475 / exptime475
    mab475 = -2.5 * np.log10(fnujy475 / 3631.)
        
    # do photometry on F814W image
    flux814 = [0.]
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data814, aperture)
        flux814.append(phot_table['aperture_sum'][0])

    photfnu814 = header814['PHOTFNU']
    exptime814 = header814['EXPTIME']
    fnujy814 = np.array(flux814) * photfnu814 / exptime814
    mab814 = -2.5 * np.log10(fnujy814 / 3631.)
        
    # do photometry on F160W image
    flux160 = [0.]
    for radius in radii:
        aperture = CircularAperture(positions, radius)
        phot_table = aperture_photometry(data160, aperture)
        flux160.append(phot_table['aperture_sum'][0])

    photfnu160 = header160['PHOTFNU']
    exptime160 = header160['EXPTIME']
    fnujy160 = np.array(flux160) * photfnu160 / exptime160
    mab160 = -2.5 * np.log10(fnujy160 / 3631.)
        
    # set up the plot
    ax = fig.add_subplot(1,1,1)

    # plot lines connecting photometric points for each aperture
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([fnujy475[i+1]-fnujy475[i], fnujy814[i+1]-fnujy814[i], fnujy160[i+1]-fnujy160[i]])
        ax.plot(wv, fx, color='0.75')

    # plot symbols for each photometric point with color scaled to aperture radius
    for i in range(0,len(radii)):
        wv = np.array([475.,814.,1600])
        fx = np.array([fnujy475[i+1]-fnujy475[i], fnujy814[i+1]-fnujy814[i], fnujy160[i+1]-fnujy160[i]])
        cax = ax.scatter(wv, fx, c=np.array([radii[i],radii[i],radii[i]]), vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2, marker='s')
        
    # set plot parameters
    cbar = fig.colorbar(cax)
    cbar.set_label('radius [pixels]', fontsize=18)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylim([3e-7,5e-5])
    plt.ylabel('flux [Jy]', fontsize=18)
    plt.xlim([300.,2000.])
    plt.xlabel(r'wavelength [$\AA$]', fontsize=18)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'ad_phot_sed.pdf')
