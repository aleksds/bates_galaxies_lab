# Joshua Rines
# 20161029
#
# the main goal of this code is to create color vs color plots comparing subtractive flux annuli in three filters
# the plot generated with this code will plot U-V color vs V-J color
    #it is our hope to use this plot for comparison to Williams (2009) to make determinations regarding the SFR processes in different regions of the galaxy
#it would also be good in the future to include code to create these plots for any galaxy. it is currently set up to only run for one galaxy, and we can change the galaxy being analyzed by changing the centroid locations.  this would be a simple update

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
import math

# define the directory that contains the images
dir = os.environ['HSTDIR']

#setting up arrays with three elements, all zeros - placeholders
wavelengths = [4,8,1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3386.5
ycen = 3503.2
dx = 100
dy = 100

# define the radii to be used for aperture photometry
radii = np.arange(40)+1

#make an array for the calculation of the area of each bagel (annulus)
area = [0 for x in range(len(radii))]

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# create a PDF file for the plots    
with PdfPages('jr_color_color.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir+'J0905_final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']

        #define positions for photometry
        positions = [(xcen, ycen)]

        #do photometry on images
        #convert to proper units
        for j in range(0,len(radii)):
            aperture = CircularAperture(positions, radii[j])
            phot_table = aperture_photometry(data[i], aperture)
            flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])
            if j == 0:
                subflux[i,j] = flux[i,j]
            else:
                subflux[i,j] = flux[i,j]-flux[i,j-1]
        
        #set up the plot
        ax = fig.add_subplot(1,1,1)
        cax = ax.scatter(-2.5*np.log10(subflux[1] / subflux[2]), -2.5*np.log10(subflux[0] / subflux[1]), c=radii, vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2,)

        #finding magnitudes for M/L ratio
        mag475 = -2.5*np.log10(flux[0,13])
        mag814 = -2.5*np.log10(flux[1,13])
        color814 = mag475-mag814

        color = -2.5*np.log10(subflux[1] / subflux[2])


    # set plot parameters
    cbar = fig.colorbar(cax)
    cbar.set_label('Radius [pixels]', fontsize=18)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #plt.ylim([10**(-1),1e1])
    plt.title('U-V color vs. V-J color, J0905', fontsize=18)
    plt.ylabel('U-V', fontsize=14)
    #plt.xlim([10**(-1),1e1])
    plt.xlabel('V-J', fontsize=14)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_color_color.pdf')
