#write code that will look into sdss files and write loop for each galaxy figure out each ra and dec
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import glob
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii

#read in flux riight before graphing(dont store just use and rewrite)

dir = os.environ['SDSSDIR']
file = glob.glob(dir+'*.fits')

#read in the relevant sdss files
# speed of light in Angstroms per second
c = 3E5

dir2 = os.environ['HIRESDIR']
# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir2+'gal_info.txt')
gal = gal_info['gal']
vcen = gal_info['vcen']
z = gal_info['zem']

#wavelengths:
w_int = ['MgII2796','MgII2803', 'MgI2852', 'FeII2600']
w_rest = np.array([2796, 2803 , 2852, 2600])
w_obs = np.zeros([len(gal), len(w_rest)])
#define w_obs:
for j in range(0, len(gal)):
    for i in range(0, len(w_rest)):
       w_obs[j][i] = w_rest[i]*(1+z[j])
    print(w_obs)

with PdfPages('boop') as pdf:
    for j in range(0, len(gal)):
        datafile = dir2+gal[j]+'/'+gal[j]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave2 = np.array([data['wv'] * u.AA])
        flux2 = data['norm']
        hdulist = fits.open(file[i])
        coeff0 = hdulist[0].header['COEFF0']
        coeff1 = hdulist[0].header['COEFF1']
        flux = hdulist[1].data
        model = hdulist[1].data['model']
        npix = len(flux)
        index = np.arange(npix)
        #wavelength = np.array([10.**(coeff0 + coeff1*index)])
        loglam = hdulist[1].data['loglam']
        wavelength = np.array([10.**loglam])
        flux_cor = np.zeros(len(flux))
        for i in range(0, len(flux)):
            flux_cor[i] = flux[i][0]/model[i]

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.plot((wave2[0]/(1+z[j])), flux2)
        plt.plot((wavelength[0]/(1+z[j])), flux_cor)
        
        plt.title("%s Galaxy" % (gal[j]))
        ax.set_ylim(0, 2)
        ax.set_xlim(2700,2900)
        pdf.savefig()
        plt.close()
        
        fig = plt.figure()
        for i in range(0, len(w_rest)):
            #fig = plt.figure()
            vel_kms = ((wavelength[0]-w_obs[j][i])/w_obs[j][i]) * c
            vel_kms2 = ((wave2[0]-w_obs[j][i])/w_obs[j][i]) * c
            ax = fig.add_subplot(4,1,1+i)
            plt.plot(vel_kms2, flux2, label='HIRES')
            plt.plot(vel_kms, flux_cor, label='SDSS')
            plt.title("%s Galaxy with Velocity of %s Spectrum" % (gal[j], w_int[i]))
            plt.ylabel('Flux')
            ax.set_xlim(-3000,500)
            ax.set_ylim(0, 2)
            plt.tight_layout()
        plt.xlabel('Velocity')
        pdf.savefig()
        plt.close()

    os.system("open %s &" % 'boop')
