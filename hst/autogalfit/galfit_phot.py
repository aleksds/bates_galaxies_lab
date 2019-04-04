# Aleks Diamond-Stanic
# 20190404
# Goal: Starting with a suite of suite of output fits files from GALFIT, perform aperture photometry of the data and model to measure total and residual flux

import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry

# the input directory that include the output fits files from GALFIT
dir = sys.argv[1]

files = glob.glob(dir+'/*output.fits')

xcen = 50
ycen = 50
radii = np.arange(40)+1
flux = np.zeros([len(files),len(radii)])
func = np.zeros([len(files),len(radii)])

name = 'galfit_phot_'+dir+'.pdf'

with PdfPages(name) as pdf:   
    for i in range(0, len(files)):
        hdu = fits.open(files[i])
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
        data, header = hdu[0].data, hdu[0].header
        dunc, unc_header = hdu[12].data, hdu[12].header
        for j in range(0,len(radii)):
            aperture = CircularAperture([xcen, ycen], radii[j])
            phot_table = aperture_photometry(data, aperture)
            flux[i,j] = phot_table['aperture_sum'][0]
            func_table = aperture_photometry(dunc, aperture)
            func[i,j] = func_table['aperture_sum'][0]

            
        fig = plt.figure()

        plt.scatter(radii, flux[i]/1e5)
        plt.errorbar(radii,flux[i]/1e5,yerr=func[i]/1e5)
        plt.xlabel('Radius [pixels]', fontsize=14)
        plt.ylabel('Flux [image units]', fontsize=14)
        plt.title(header['LOGFILE'])
        
        pdf.savefig()
        plt.close()

os.system('open %s &' % name)

