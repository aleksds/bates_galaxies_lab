# Imports
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from matplotlib.backends.backend_pdf import PdfPages

# This is a bitmask handling object from the DAP source code
from mangadap.dapcube import DAPCubeBitMask

dapdir = os.environ['DAPDIR']


# define plate and ifu
redshift = [0.0284, 0.0284, 0.0284, 0.0284, 0.0284, 0.0284]#[0.01874]
plate = [7495, 7495, 7495, 7495, 7495, 7495]#[7443]
ifu = [12702, 12702, 12702, 12702, 12702, 12702]#[12704]
a=[38, 27, 51, 50, 32, 27]
b=[37, 29, 27, 43, 21, 59]

for i in range(0, len(plate)):

    # find the fits files
    hdu_maps_file = glob.glob(dapdir+'SPX-GAU-MILESHC/'+str(plate[i])+'/'+str(ifu[i])+'/'+'manga-'+str(plate[i])+'-'+str(ifu[i])+'-MAPS-SPX-GAU*')
    hdu_cube_file = glob.glob(dapdir+'SPX-GAU-MILESHC/'+str(plate[i])+'/'+str(ifu[i])+'/'+'manga-'+str(plate[i])+'-'+str(ifu[i])+'-LOGCUBE-SPX-GAU*')
    
    # Open the fits file
    hdu_maps = fits.open(hdu_maps_file[0])
    hdu_cube = fits.open(hdu_cube_file[0])
    
    
    #hdu_maps = fits.open(dir+'manga-7443-12704-MAPS-SPX-GAU-MILESHC.fits.gz')
    #hdu_cube = fits.open(dir+'manga-7443-12704-LOGCUBE-SPX-GAU-MILESHC.fits.gz')
    
    # Get the S/N per bin from the MAPS file
    snr = np.ma.MaskedArray(hdu_maps['BIN_SNR'].data, mask=hdu_maps['BINID'].data < 0)
    
    # Select the bin/spaxel with the highest S/N
    k = np.ma.argmax(snr.ravel())
    n = hdu_maps['BIN_SNR'].data.shape[0] # Number of pixels in X and Y

    # Get the pixel coordinate
    l=a[i]
    j=b[i]
    
    # Declare the bitmask object to mask selected pixels
    bm = DAPCubeBitMask()
    wave = hdu_cube['WAVE'].data
    flux = np.ma.MaskedArray(hdu_cube['FLUX'].data[:,j,l],
                                mask=bm.flagged(hdu_cube['MASK'].data[:,j,l],
    				[ 'IGNORED', 'FLUXINVALID', 'IVARINVALID', 'ARTIFACT' ]))
    model = np.ma.MaskedArray(hdu_cube['MODEL'].data[:,j,l],
                                 mask=bm.flagged(hdu_cube['MASK'].data[:,j,l], 'FITIGNORED'))
    stellarcontinuum = np.ma.MaskedArray(
                            hdu_cube['MODEL'].data[:,j,i] - hdu_cube['EMLINE'].data[:,j,l]
                                - hdu_cube['EMLINE_BASE'].data[:,j,l],
                                 mask=bm.flagged(hdu_cube['MASK'].data[:,j,l], 'FITIGNORED'))
    emlines = np.ma.MaskedArray(hdu_cube['EMLINE'].data[:,j,l],
                                   mask=bm.flagged(hdu_cube['EMLINE_MASK'].data[:,j,l],
                                                   'ELIGNORED'))
    resid = flux-model-0.5
    
    # create a PDF file for the plots    
    with PdfPages('dap_maps_logcube_'+str(plate[i])+'_'+str(ifu[i])+'_'+str(a[i])+'_'+str(b[i])+'.pdf') as pdf:
    
        fig = plt.figure()
    
        ax = fig.add_subplot(3,1,1)
        model_median= np.median(model)
        ax.step(wave, flux, where='mid', color='k', lw=0.5)
        ax.plot(wave, model, color='r', lw=1)
        #ax.plot(wave, stellarcontinuum, color='g', lw=1)
        #ax.plot(wave, emlines, color='b', lw=1)
        ax.step(wave, resid, where='mid', color='0.5', lw=0.5)
        plt.xlim([3600,10500])
        plt.ylim([-0.1*model_median,3*model_median])
        plt.xlabel('Wavelength [$\AA$]',fontsize=10)
        plt.ylabel('$f_{\lambda}$ [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',fontsize=8)
        ax.text(6000,0,str(a[i])+','+str(b[i]))
        pdf.savefig()
        plt.close()
        fig=plt.figure()
        
        lambda_zero=3727*(1.+redshift[i])
        lmin = lambda_zero-30
        lmax = lambda_zero+30
        ymax = np.max(model[(wave > lmin) & (wave < lmax)])
        
        ax = fig.add_subplot(3,1,1)
        ax.step(wave, flux, where='mid', color='k', lw=0.5)
        ax.plot(wave, model, color='r', lw=1)
        #ax.plot(wave, stellarcontinuum, color='g', lw=1)
        #ax.plot(wave, emlines, color='b', lw=1)
        plt.xlim([lambda_zero-30,lambda_zero+30])
        plt.ylim([-.05,ymax+0.1])
        plt.xlabel('Wavelength [$\AA$]',fontsize=10)
        plt.ylabel('$f_{\lambda}$ [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',fontsize=8)
        ax.text(lambda_zero,0,str(a[i])+','+str(b[i]))
        
        
        ax = fig.add_subplot(3,1,3)
        vel=3E5*(wave-lambda_zero)/lambda_zero
        ax.step(vel, flux, where='mid', color='k', lw=0.5)
        ax.plot(vel, model, color='r', lw=1)
        #ax.plot(vel, stellarcontinuum, color='g', lw=1)
        #ax.plot(vel, emlines, color='b', lw=1)
        plt.xlim([-1000,1000])
        plt.ylim([-.05,ymax+0.1])
        plt.xlabel('Velocity (km s$^{-1}$)',fontsize=10)
        plt.ylabel('$f_{\lambda}$ [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',fontsize=8)
        ax.text(0,0,str(a[i])+','+str(b[i]))
        
        pdf.savefig()
        plt.close()

    name = 'dap_maps_logcube_'+str(plate[i])+'_'+str(ifu[i])+'_'+str(a[i])+'_'+str(b[i])+'.pdf'
    os.system('open %s &' % name)
