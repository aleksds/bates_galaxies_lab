# Joshua Rines
# 20161029
#
# the main goal of this code is to perform photometry in different
# apertures in three different images for a galaxy and then plots flux
# vs wavelegnth (i.e., a spectral energy distribution) for each aperture

#includes conversion of flux from image units to Jy
#includes dividing flux by area

#would be good to make this compatible with any galaxy instead of needing to change the centroid values for each galaxy

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
xcen = 3572.9
ycen = 3339.1
dx = 100
dy = 100

# define the radii to be used for aperture photometry
radii = np.array([1,2,3,4,5,8,10,15,20,25,30,35,40])

#make an array for the calculation of the area of each bagel (annulus)
area = [0 for x in range(len(radii))]

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# create a PDF file for the plots    
with PdfPages('jr_phot_sed.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir+'J1107_final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']

        #define positions for photometry
        positions = [(xcen, ycen)]

        #do photometry on images
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

    # plot lines connecting photometric points for each aperture
    for j in range(0,len(radii)-1):
        wv = np.array([475.,814.,1600.])
        fx = (subflux[0:len(wv),j])/(area[j])
        ax.plot(wv, fx, color='0.75')
        
        # plot symbols for each photometric point with color scaled to aperture radius
        cax = ax.scatter(wv, fx, c=np.array([radii[j],radii[j],radii[j]]), vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2, marker='s')


    # set plot parameters
    cbar = fig.colorbar(cax)
    cbar.set_label('Radius [pixels]', fontsize=14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.title('Subflux vs. Wavelength, J1107', fontsize=18)
    #plt.ylim([1e4,6e5])
    plt.ylabel('Flux/Area (Jy/pixel)', fontsize=14)
    plt.xlim([300.,2000.])
    plt.xlabel(r'Wavelength [$\AA$]', fontsize=14)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_phot_sed.pdf')
