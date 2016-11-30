# Josh Rines
# 20161014

# the preliminary goals of this code include the following:

# (1) read in three images taken in three different filters
# (2) do aperature photometry on each image
# (3) use measurements in three different aperatures to make a plot of flux vs wavelength for each aperature

# we have since then gone on to make three plots with this code; all of which plot the the data from three different filters showing F475W as blue, F814W as green, and F160W as red, as listed below:
# (1) log of subtractive flux vs wavelength
# (2) subtractive flux vs aperature radius
# (3) flux vs wavelength


# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

# define a function to plot "postage stamp" images
def plot_image(stamp):
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)

# define the directory that contains the images
dir = os.environ['HSTDIR']

# specify the position of the science target and the size of the
# region around the science target to consider
xcen = 3386.5
ycen = 3503.2
dx = 100
dy = 100

# create a PDF file for the plots    
with PdfPages('jr_flux_plots.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W', 'F814W', 'F160W']
    filter = [475, 814, 1600]
    colors = ['bo', 'go', 'ro']
    colorz = ['b','g','r']
    
    for i in range (0, len(collection)):

        # read in the images
        file = glob.glob(dir+'J0905_final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data, header = hdu[0].data, hdu[0].header
    
        # create a "postage stamp" image centered on the science target
        stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]
    
        # plot the "postage stamp"
        # ax = fig.add_subplot(2,3,1)
        # plot_image(stamp475)
        # plt.suptitle(header475['TARGNAME'])
        # plt.tick_params(axis='both', which='major', labelsize=8)
    
        # set up for plots, establishing "fence post" correction and defining list lengths
        positions = [(xcen, ycen)]
        radii = np.arange(dx)+1
        wavelength = np.arange(dx)+1
        subflux = np.arange(dx)+1
        flux = []
        for radius in radii:
            aperture = CircularAperture(positions, radius)
            phot_table = aperture_photometry(data, aperture)
            flux.append(phot_table['aperture_sum'][0])
            #print (phot_table)
        for j in range (0, len(radii)):
            wavelength[j] = filter[i]
            if j == 0:
                subflux[j] = flux[0]
            else:
                subflux[j] = flux[j]-flux[j-1]
        
        labeling = ['F475W','814W','160W']
            
        # plot log of subflux vs wavelength for three filters  
        ax = fig.add_subplot(3,1,1)
        ax.plot(wavelength, np.log10(subflux), colors[i])
        plt.title('Log Subflux vs. Wavelength (three filters), J0905')
        plt.xlim([400,1700])
        plt.ylim([3,6])
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('log Subflux (dv)')
        plt.tick_params(axis='both', which='major', labelsize=8)
        fig.tight_layout()

        # plot the subtractive flux vs aperature radius
        bx = fig.add_subplot(3,1,2)
        bx.plot(radii, subflux, colorz[i])
        plt.title('Subtractive Flux vs. Aperture Radius, J0905')
        plt.xlabel('Aperature Radius (pixels)')
        plt.ylabel('Subflux (dv)')
        plt.tick_params(axis='both', which='major', labelsize=8)
        fig.tight_layout()

        # plot the flux vs radius
        cx = fig.add_subplot(3,1,3)
        cx.plot(radii, flux, colorz[i])
        plt.title('Flux vs. Radius, J0905')
        plt.xlabel('Aperature Radius (pixels)')
        plt.ylabel('Flux (dv)')
        plt.tick_params(axis='both', which='major', labelsize=8)
        fig.tight_layout()


    
    pdf.savefig()
    plt.close()


    os.system('open %s &' % 'jr_flux_plots.pdf')
