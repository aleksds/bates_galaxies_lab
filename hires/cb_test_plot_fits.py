from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


#read relevant quantities
def get_quantities(filepath):
    #open file given filepath
    fits_file = fits.open(filepath)

    flux_values = fits_file['COADD'].data['flux']
    model_val = fits_file['COADD'].data['model']
    HBWAV = fits_file['SPZLINE'].data['LINEWAVE'][15]*(1+fits_file['SPZLINE'].data['LINEZ'][15])

    minw = fits_file['SPECOBJ'].data['WAVEMIN']
    maxw = fits_file['SPECOBJ'].data['WAVEMAX']

    data_n = len(flux_values)

    wave = 10**fits_file['COADD'].data['loglam'] #wavelegnths corresponding to each flux value

    return [wave,flux_values,model_val,minw,maxw,HBWAV]

#plot
def pdfplot_flux_deffit(data,outputfile):
    with PdfPages(outputfile) as pdf:
        wave,flux_values,model_val,minw,maxw,HBWAV = data[0],data[1],data[2],data[3],data[4],data[5]

        fig, ax = plt.subplots(dpi=1200)

        plt.axvline(x=HBWAV, color='blue', label='Hb', alpha=0.4)

        ax.plot(wave,flux_values,color='black',linewidth=0.3)
        ax.plot(wave,model_val,color='red',linewidth=0.3, alpha=0.7)

        plt.grid(True)
        plt.xlim(minw,maxw)
        ax.xaxis.set_major_locator(plt.MultipleLocator(500))
        plt.xlabel("$\AA ngstr \ddot{o} ms$")
        plt.ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")
        plt.legend()

        pdf.savefig()
        plt.close('all')

#vars
filepath = "/Users/cbradna/Documents/spec-0761-54524-0409.fits"
outputfile = "plotpdftest.pdf"

#Run routine
data = get_quantities("/Users/cbradna/Documents/spec-0761-54524-0409.fits")
pdfplot_flux_deffit(data,outputfile)