#code to read galfit output files and extract the desired information (integrated magnitudes)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']
number = ['01','02']
mags = np.zeros([2,12,2]) #order: psf models then sersic models, normal order for galaxies and filters

for m in range(0,2):
    for w in range(0,12):
        for i in range(0,2):

            file = galaxies[w]+'_'+filters[i]+'_coarse.galfit.'+number[m]+'.band'
            with open(file) as f:
                content = f.readlines()
            mags[m][w][i] = np.float(content[47][4:10])

psfyvals = np.zeros(12)
psfxvals = np.zeros(12)
sersicyvals = np.zeros(12)
sersicxvals = np.zeros(12)
for w in range(0,12):
    psfxvals[w] = mags[0][w][1]
    sersicxvals[w] = mags[1][w][1]
    psfyvals[w] = mags[0][w][0] - mags[0][w][1]
    sersicyvals[w] = mags[1][w][0] - mags[1][w][1]
    
with PdfPages('psf&sersic_color_mag.pdf') as pdf:   
    plt.figure()
    plt.scatter(psfxvals,psfyvals, label='PSF')
    plt.scatter(sersicxvals,sersicyvals, label='SERSIC')
    plt.xlabel('F814W')
    plt.ylabel('F475W - F814W')
    plt.title('Color Magnitude Plot (with psf and sersic fits)')
    plt.legend(loc='upper right')
    pdf.savefig()
    plt.close()
