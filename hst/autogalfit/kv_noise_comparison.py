#This is going to be the code for the image comparison between galfit sigma images and our created sigma images
#the purpose of this piece of code is to find out if the uncertainty calculations are right by doing a comparison
#first recreating our sigma images
#taking code from ad_monster.py

#relevant packages

import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from scipy import integrate
from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
import img_scale
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry


#now reading and storing galfit images
#essentially stealing code from ad_galfit_phot.py & ad_monster.py

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

pixNoise = np.zeros([len(filters), width, width])

# the code requires the user to specify an input directory that has output fits files from GALFIT
dir = sys.argv[1]
gal = sys.argv[2]

# defining galaxies
# galaxies =['j0826', 'j0901','j0905', 'j0944', 'j1107', 'j1219', 'j1341', 'j1506', 'j1558', 'j1613', 'j2116','j2140']

# find all of the relevant fits files
files = glob.glob(dir + '/' + gal + '*output.fits')

# define relevant variables
xcen = 50
ycen = 50
fmin = -2
fmax = 12
# radii = np.arange(40)+1
# radii = 10.**(np.arange(21)/20*1.9-0.3)
# radii = np.array([1.,2.,3.,4.,5.,10.,15.,20.,25.,30.,35.,40.,45.])
radii = np.array([5., 10., 20., 30., 40.])

# arrays for v band
vdat_flux = np.zeros([len(files), len(radii)])
vunc_flux = np.zeros([len(files), len(radii)])
vmod_flux = np.zeros([len(files), len(radii)])
vres_flux = np.zeros([len(files), len(radii)])
vres_sub = np.zeros([len(files), len(radii)])
vdat_sub = np.zeros([len(files), len(radii)])

# arrays for u band
udat_flux = np.zeros([len(files), len(radii)])
uunc_flux = np.zeros([len(files), len(radii)])
umod_flux = np.zeros([len(files), len(radii)])
ures_flux = np.zeros([len(files), len(radii)])
ures_sub = np.zeros([len(files), len(radii)])
udat_sub = np.zeros([len(files), len(radii)])

# arrays for j band
jdat_flux = np.zeros([len(files), len(radii)])
junc_flux = np.zeros([len(files), len(radii)])
jmod_flux = np.zeros([len(files), len(radii)])
jres_flux = np.zeros([len(files), len(radii)])
jres_sub = np.zeros([len(files), len(radii)])
jdat_sub = np.zeros([len(files), len(radii)])

# let's produce plots that will be saved in a pdf file
name = 'galfit_phot_' + dir + '_' + gal + '.pdf'

with PdfPages(name) as pdf:
    # loop over each output fits file from GALFIT
    for i in range(0, len(files)):

        # read in the image
        hdu = fits.open(files[i])

        # these images have 18 different extensions:
        # hdu[0] is INPUT_V
        # hdu[1] is INPUT_U
        # hdu[2] is INPUT_J
        # hdu[3] is MODEL_V
        # hdu[4] is MODEL_U
        # hdu[5] is MODEL_J
        # hdu[6] is RESIDUAL_V
        # hdu[7] is RESIDUAL_U
        # hdu[8] is RESIDUAL_J
        # hdu[9] is COMPONENT_1_sersic_V
        # hdu[10] is COMPONENT_1_sersic_U
        # hdu[11] is COMPONENT_1_sersic_J
        # hdu[12] is SIGMA_V
        # hdu[13] is SIMGA_U
        # hdu[14] is SIGMA_J
        # hdu[15] is PSF_V
        # hdu[16] is PSF_U
        # hdu[17] is PSF_J

        # F814W corresponds roughly to rest-frame V-band
        vdat, vdat_head = hdu[0].data, hdu[0].header
        vunc, vunc_head = hdu[12].data, hdu[12].header
        vmod, vmod_head = hdu[3].data, hdu[3].header
        vres, vres_head = hdu[6].data, hdu[6].header

        # F475W corresponds roughly to rest-frame U-band
        udat, udat_head = hdu[1].data, hdu[1].header
        uunc, uunc_head = hdu[13].data, hdu[13].header
        umod, umod_head = hdu[4].data, hdu[4].header
        ures, ures_head = hdu[7].data, hdu[7].header

        # F160W corresponds roughly to rest-frame J-band
        jdat, jdat_head = hdu[2].data, hdu[2].header
        junc, junc_head = hdu[14].data, hdu[14].header
        jmod, jmod_head = hdu[5].data, hdu[5].header
        jres, jres_head = hdu[8].data, hdu[8].header

        # loop over each aperture
        for j in range(0, len(radii)):

            # define the circular aperture
            aperture = CircularAperture([xcen, ycen], radii[j])

            # perform aperture photmetry in that aperture and save the flux value
            vdat_table = aperture_photometry(vdat, aperture)
            vdat_flux[i, j] = vdat_table['aperture_sum'][0]
            udat_table = aperture_photometry(udat, aperture)
            udat_flux[i, j] = udat_table['aperture_sum'][0]
            jdat_table = aperture_photometry(jdat, aperture)
            jdat_flux[i, j] = jdat_table['aperture_sum'][0]

            # repeat for the sigma image
            vunc_table = aperture_photometry(vunc, aperture)
            vunc_flux[i, j] = vunc_table['aperture_sum'][0]
            uunc_table = aperture_photometry(uunc, aperture)
            uunc_flux[i, j] = uunc_table['aperture_sum'][0]
            junc_table = aperture_photometry(junc, aperture)
            junc_flux[i, j] = junc_table['aperture_sum'][0]

            # repeat for the model image
            vmod_table = aperture_photometry(vmod, aperture)
            vmod_flux[i, j] = vmod_table['aperture_sum'][0]
            umod_table = aperture_photometry(umod, aperture)
            umod_flux[i, j] = umod_table['aperture_sum'][0]
            jmod_table = aperture_photometry(jmod, aperture)
            jmod_flux[i, j] = jmod_table['aperture_sum'][0]

            # repeat for the residual image
            vres_table = aperture_photometry(vres, aperture)
            vres_flux[i, j] = vres_table['aperture_sum'][0]
            ures_table = aperture_photometry(ures, aperture)
            ures_flux[i, j] = ures_table['aperture_sum'][0]
            jres_table = aperture_photometry(jres, aperture)
            jres_flux[i, j] = jres_table['aperture_sum'][0]


        # the following figures produce the galfit sigma images in a logscale
        fig = plt.figure()
        logpiximg475 = img_scale.log(uunc, scale_min=1, scale_max=500)
        plt.imshow(logpiximg475)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F475W (Log)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        logpiximg814 = img_scale.log(vunc, scale_min=1, scale_max=500)
        plt.imshow(logpiximg814)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F814W (Log)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        logpiximg160 = img_scale.log(junc, scale_min=1, scale_max=500)
        plt.imshow(logpiximg160)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F160W (Log)')
        pdf.savefig()
        plt.close()

        #these produce linear scale galfit images

        fig = plt.figure()
        linpiximg475 = img_scale.linear(uunc, scale_min=1, scale_max=500)
        plt.imshow(linpiximg475)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F475W (Linear)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        linpiximg814 = img_scale.linear(vunc, scale_min=1, scale_max=500)
        plt.imshow(linpiximg814)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F814W (Linear)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        linpiximg160 = img_scale.linear(junc, scale_min=1, scale_max=500)
        plt.imshow(linpiximg160)
        plt.colorbar()
        plt.title('J2140_' + 'pixnoise_' + 'F160W (Linear)')
        pdf.savefig()
        plt.close()

        #testing the images produced by our own creation here in linear
        fig = plt.figure()
        noiseimg475 = img_scale.linear(pixNoise[0], scale_min=1, scale_max=500)
        plt.imshow(noiseimg475)
        plt.colorbar()
        plt.title(galaxies[w].name + 'pixNoise' + 'F475W')
        pdf.savefig()
        plt.close()

        # now coding the comparison plot between galfit sigma images and pixnoise for filter F475W
        ratio475 = pixNoise[0]/uunc
        fig = plt.figure()
        ratioimg475 = img_scale.linear(ratio475, scale_min=1, scale_max=500)
        plt.imshow(ratioimg475)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F475W (Linear)')
        pdf.savefig()
        plt.close()

        #this will be the same plot but now in log form
        fig = plt.figure()
        logratioimg475 = img_scale.log(ratio475, scale_min=1, scale_max=500)
        plt.imshow(logratioimg475)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F475W (Log)')
        pdf.savefig()
        plt.close()

        #filter F814W
        ratio814 = pixNoise[1] / vunc
        fig = plt.figure()
        ratioimg814 = img_scale.linear(ratio814, scale_min=1, scale_max=500)
        plt.imshow(ratioimg814)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F814W (Linear)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        logratioimg814 = img_scale.log(ratio814, scale_min=1, scale_max=500)
        plt.imshow(logratioimg814)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F814W (Log)')
        pdf.savefig()
        plt.close()

        # filter F160W
        ratio160 = pixNoise[2] / junc
        fig = plt.figure()
        ratioimg160 = img_scale.linear(ratio160, scale_min=1, scale_max=500)
        plt.imshow(ratioimg160)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F160W (Linear)')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        logratioimg160 = img_scale.log(ratio160, scale_min=1, scale_max=500)
        plt.imshow(logratioimg160)
        plt.colorbar()
        plt.title('J2140 Pixnoise Galfit Sigma Ratio for F160W (Log)')
        pdf.savefig()
        plt.close()


        size = len(ratioimg475)
        xpos = np.zeros([size,size])
        ypos = np.zeros([size,size])

        for j in range(0, size):
            for k in range(0, size):
                xpos[j, k] = k
                ypos[j, k] = j

        xpix = xpos - np.median(xpos)
        ypix = ypos - np.median(ypos)

        fig = plt.figure()
        radius = np.sqrt((xpix-xcen)**2+(ypix-ycen)**2)
        ratiopix475 = ratioimg475.flatten()
        plt.scatter(radius,ratiopix475)
        plt.xlim(0,73)
        plt.ylim(0,0.01)
        plt.title('Ratio Pixel Values vs Radius for F475W (Linear)')
        plt.ylabel('Ratio Pixel Values')
        plt.xlabel('Radius')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        radius = np.sqrt((xpix - xcen) ** 2 + (ypix - ycen) ** 2)
        logratiopix475 = logratioimg475.flatten()
        plt.scatter(radius, logratiopix475)
        #plt.xlim(0, 73)
        #plt.ylim(0, 0.01)
        plt.title('Ratio Pixel Values vs Radius for F475W (Log)')
        plt.ylabel('Ratio Pixel Values')
        plt.xlabel('Radius')
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        radius = np.sqrt((xpix - xcen) ** 2 + (ypix - ycen) ** 2)
        ratiopix475 = ratioimg475.flatten()
        plt.scatter(radius, np.log(ratiopix475))
        #plt.xlim(0, 73)
        #plt.ylim(0, 0.01)
        plt.title('Ratio Pixel Values vs Radius for F475W (Log2)')
        plt.ylabel('Ratio Pixel Values')
        plt.xlabel('Radius')
        pdf.savefig()
        plt.close()

# open the pdf file
os.system('open %s &' % name)