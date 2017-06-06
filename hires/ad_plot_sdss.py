# Aleks Diamond-Stanic, 20170606
# code to read in and plot an SDSS spectrum

import os
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

dir = os.environ['SDSSDIR']
file = dir + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)

plt.plot(hdulist[5].data)

filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(hdulist[5].data)

    pdf.savefig()
    plt.close()
    
    os.system("open %s &" % filename)
