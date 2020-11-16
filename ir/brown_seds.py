import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter

projpath = os.getcwd()+'/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/'

files = glob.glob(projpath+'*.dat')


name = 'brown_seds.pdf'

with PdfPages(name) as pdf:
    for i in range(0,len(files)):
        fig = plt.figure()

        ax = fig.add_subplot(111)

        table = ascii.read(files[i])

        wave = np.array(table['col1'])
        flam = np.array(table['col2'])
        flux = wave * flam

        plt.plot(wave/1e4, flux)

        loc = files[0].find('spec')

        gal = files[i][loc-9:loc-1]
        print(gal)
        
        plt.title(gal)

        prange = wave > 3e4
        fmin = np.min(flux[prange])
        fmax = np.max(flux[prange])
        
        plt.yscale('log')
        plt.ylim([fmin/2,fmax*2])

        plt.xscale('log')
        #plt.xlim([3e4,4e5])
        plt.xlim(3,40)

        plt.xlabel(r'Wavelength [$\mu$m]')
        plt.ylabel(r'Flux $\lambda f_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$]')

        ax.set_xticks([3, 5, 8, 10, 20, 30])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        
        labels = [item.get_text() for item in ax.get_xticklabels()]
        #labels[1] = 'Testing'
        labels[0] = '3'
        labels[1] = '5'
        labels[2] = '8'
        labels[3] = '10'
        labels[4] = '20'
        labels[5] = '30'
        #print(labels)
        


        pdf.savefig()
        plt.close()

os.system('open %s &' % name)
