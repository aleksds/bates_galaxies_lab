# Aleks Diamond-Stanic, 20170606
# code to read in and plot an SDSS spectrum

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

dir = os.environ['SDSSDIR']
file = dir + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)

coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']

flux = hdulist[1].data

npix = len(flux)
index = np.arange(npix)
wavelength = 10.**(coeff0 + coeff1*index)

filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(wavelength, flux)
    ax.set_ylim(0., 10.)
    pdf.savefig()
    plt.close()
    
    os.system("open %s &" % filename)
