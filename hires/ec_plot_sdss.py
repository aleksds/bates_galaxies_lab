#Eve Cinquino, 20170606
#Read in an SDSS spectrum and make a plot of flux versus wavelength
#
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator

dir = os.environ['SDSSDIR']
file = dir + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)

coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']

flux = hdulist[1].data

npix = len(flux)
index = np.arange(npix)
wavelength = 10.**(coeff0 + coeff1*index)

z = 0.60329
w_int = ['MgII2796','MgII2803', 'MgI2852', 'FeII2600']
w_rest = np.array([2796, 2803, 2852, 2600])
w_obs = []
for i in range(0, len(w_rest)):
    temp = w_rest[i]*(1+z)
    w_obs.append(temp)
print(w_obs)

minorLocator = AutoMinorLocator()
filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:
    #prints the first graph 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(wavelength, flux)
    for i in range(0, len(w_rest)):
        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    ax.set_ylim(0., 15.)
    ax.set_xlim(3800,9000)
    pdf.savefig()
    plt.close()

    #prints each of the zoomed in graphs
    
    for i in range(0, len(w_int)):
        # read in the spectrum
        fig2 = plt.figure()
        ax = fig2.add_subplot(1,1,1)
        ax.plot([(1, 2), (3, 4)], [(4, 3), (2, 3)])
        plt.title(w_int[i])
        ax.set_xlim(w_obs[i]-200, w_obs[i]+200)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.plot(wavelength, flux)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+10., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-10., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        pdf.savefig()
        plt.close()
    
    os.system("open %s &" % filename)
    
