# Aleks Diamond-Stanic 20170706
# based on code written by John Moustakas

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import fitsio
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM

from matplotlib.backends.backend_pdf import PdfPages

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Read the original photometry, the fitting results, and the K-corrections.
massdir = os.getenv('ISEDFITDIR')
photfile = 'sg_fluxtable_nm.txt'
isedfile = os.path.join(massdir, 'massprofiles_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits.gz')
kcorrfile = os.path.join(massdir, 'massprofiles_fsps_v2.4_miles_chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz')

print('Reading {}'.format(photfile))
phot = ascii.read(photfile)
phot[:2]

print('Reading {}'.format(isedfile))
ised = Table(fitsio.read(isedfile, ext=1, upper=True))
ised[:2]

print('Reading {}'.format(kcorrfile))
kcorr = Table(fitsio.read(kcorrfile, ext=1, upper=True))
kcorr[:2]

galaxy = [gg[:5] for gg in phot['ID'].data]
galaxy = np.unique(galaxy)
ngal = len(galaxy)

# Plot the individual stellar mass profiles.

nrad = 40
radpix = np.linspace(1.0, 40.0, nrad) # [pixels]
radarcsec = radpix * 0.05             # [arcsec]

mstar = ised['MSTAR_AVG'].data.reshape(ngal, nrad)
mstar_err = ised['MSTAR_ERR'].data.reshape(ngal, nrad)
redshift = phot['z'].data.reshape(ngal, nrad)[:, 0]

area = np.pi * np.insert(np.diff(radarcsec**2), 0, radarcsec[0]**2) # aperture annulus [arcsec2]

sigma = np.zeros_like(mstar)  # surface mass density [Mstar/kpc2]
radkpc = np.zeros_like(mstar) # radius [comoving kpc]

for igal in range(ngal):
    arcsec2kpc = cosmo.arcsec_per_kpc_comoving(redshift[igal]).value
    radkpc[igal, :] = radarcsec / arcsec2kpc
    areakpc2 = area / arcsec2kpc**2
    sigma[igal, :] = np.log10( 10**mstar[igal, :] / areakpc2 )

    
massrange = (8, 10.2)
sigmarange = (6, 9.6)

filename='ad_massprofiles.pdf'

with PdfPages(filename) as pdf:
    

    fig, ax = plt.subplots(3, 4, figsize=(14, 8), sharey=True, sharex=True)
    for ii, thisax in enumerate(ax.flat):
        thisax.errorbar(radarcsec, mstar[ii, :], yerr=mstar_err[ii, :], 
                        label=galaxy[ii])
        thisax.set_ylim(massrange)
        #thisax.legend(loc='upper right', frameon=False)  
        thisax.annotate(galaxy[ii], xy=(0.9, 0.9), xycoords='axes fraction', 
                        size=16, ha='right', va='top')
    fig.text(0.0, 0.5, r'$\log_{10}\, (M / M_{\odot})$', ha='center', 
             va='center', rotation='vertical')
    fig.text(0.5, 0.0, 'Radius (arcsec)', ha='center', 
             va='center')
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.tight_layout()
    
    pdf.savefig()
    plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    for igal in range(ngal):
        ax.plot(radkpc[igal, :], np.log10(np.cumsum(10**mstar[igal, :])), label=galaxy[igal])
    ax.legend(loc='lower right', fontsize=16, ncol=3, frameon=False)
    ax.set_xlabel(r'Galactocentric Radius $r_{kpc}$ (Comoving kpc)')
    ax.set_ylabel(r'$\log_{10}\, M(<r_{kpc})\ (M_{\odot})$')

    pdf.savefig()

    plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    for igal in range(ngal):
        ax.plot(radkpc[igal, :], sigma[igal, :], label=galaxy[igal])
    ax.legend(loc='upper right', fontsize=16, ncol=3, frameon=False)
    ax.set_xlabel(r'Galactocentric Radius $r_{kpc}$ (Comoving kpc)')
    ax.set_ylabel(r'$\log_{10}\, \Sigma\ (M_{\odot}\ /\ {\rm kpc}^2)$')
    ax.set_ylim(sigmarange)

    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)
