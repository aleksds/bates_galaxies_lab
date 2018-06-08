import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

mpl7_dir = os.environ['MANGADIR_MPL7'] 
drp = fits.open(mpl7_dir+'drpall-v2_4_3.fits')
drpdata = drp[1].data 


ba = drpdata.field('nsa_sersic_ba')
stellar_mass = drpdata.field('nsa_sersic_mass')
absmag = drpdata.field('nsa_elpetro_absmag')
gr= absmag[:,3] - absmag[:,4]
filename = 'histograms.pdf'


with PdfPages(filename) as pdf:

    #plot 1: b/a ratio
    fig = plt.figure()
    good = ba > 0
    print(np.count_nonzero(good))
    print(np.shape(drpdata))
    plt.hist(ba[good], color='teal', bins=np.arange(0.1,1.02,.02))
    plt.xlabel('b/a ratio')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of b/a Values in MaNGA sample')
    plt.axvline(x=0.3, color='k', linestyle='dashed')

    pdf.savefig()
    plt.close

    #plot 2: stellar mass
    fig = plt.figure()
    print(np.shape(drpdata))
    plt.hist(np.log10(stellar_mass[good]),color='teal', bins=np.arange(8.0,12.0,.2))
    plt.xlim(8,12)
    plt.xlabel('Stellar Mass (logM☉)')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of Stellar Mass in MaNGA sample')

    pdf.savefig()
    plt.close

    #plot 3: g-r
    fig = plt.figure()
    print(np.shape(drpdata))
    plt.hist(gr[good], bins=np.arange(0.0,1.0,0.05), color='teal')
    plt.xlabel('g-r')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of g-r in MaNGA sample')

    pdf.savefig()
    plt.close


    
os.system("open %s &" % filename)
    
