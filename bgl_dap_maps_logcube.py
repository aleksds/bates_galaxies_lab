# Imports
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages

# This is a bitmask handling object from the DAP source code
from mangadap.dapcube import DAPCubeBitMask

dir = '/Users/adiamond/data/MPL-5/dap/SPX-GAU-MILESHC/7977/12704/'

# Open the fits file
hdu_maps = fits.open(dir+'manga-7977-12704-MAPS-SPX-GAU-MILESHC.fits.gz')
hdu_cube = fits.open(dir+'manga-7977-12704-LOGCUBE-SPX-GAU-MILESHC.fits.gz')

# Get the S/N per bin from the MAPS file
snr = np.ma.MaskedArray(hdu_maps['BIN_SNR'].data, mask=hdu_maps['BINID'].data < 0)

# Select the bin/spaxel with the highest S/N
k = np.ma.argmax(snr.ravel())
n = hdu_maps['BIN_SNR'].data.shape[0] # Number of pixels in X and Y
# Get the pixel coordinate
j = k//n
i = k - j*n

# Declare the bitmask object to mask selected pixels
bm = DAPCubeBitMask()
wave = hdu_cube['WAVE'].data
flux = np.ma.MaskedArray(hdu_cube['FLUX'].data[:,j,i],
                            mask=bm.flagged(hdu_cube['MASK'].data[:,j,i],
				[ 'IGNORED', 'FLUXINVALID', 'IVARINVALID', 'ARTIFACT' ]))
model = np.ma.MaskedArray(hdu_cube['MODEL'].data[:,j,i],
                             mask=bm.flagged(hdu_cube['MASK'].data[:,j,i], 'FITIGNORED'))
stellarcontinuum = np.ma.MaskedArray(
                        hdu_cube['MODEL'].data[:,j,i] - hdu_cube['EMLINE'].data[:,j,i]
                            - hdu_cube['EMLINE_BASE'].data[:,j,i],
                             mask=bm.flagged(hdu_cube['MASK'].data[:,j,i], 'FITIGNORED'))
emlines = np.ma.MaskedArray(hdu_cube['EMLINE'].data[:,j,i],
                               mask=bm.flagged(hdu_cube['EMLINE_MASK'].data[:,j,i],
                                               'ELIGNORED'))
resid = flux-model-0.5

# create a PDF file for the plots    
with PdfPages('dap_maps_logcube.pdf') as pdf:

    fig = plt.figure()

    ax = fig.add_subplot(1,1,1)
    
    ax.step(wave, flux, where='mid', color='k', lw=0.5)
    ax.plot(wave, model, color='r', lw=1)
    ax.plot(wave, stellarcontinuum, color='g', lw=1)
    ax.plot(wave, emlines, color='b', lw=1)
    ax.step(wave, resid, where='mid', color='0.5', lw=0.5)

    pdf.savefig()
    plt.close()

os.system('open %s &' % 'dap_maps_logcube.pdf')
    
