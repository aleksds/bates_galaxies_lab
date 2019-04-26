from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


#open data
filepath = "/Users/cbradna/Documents/spec-0761-54524-0409.fits"
fits_file = fits.open(filepath)

#read relevant quantities
flux_values = fits_file['COADD'].data['flux']
model_val = fits_file['COADD'].data['model']
HBWAV = fits_file['SPZLINE'].data['LINEWAVE'][15]*(1+fits_file['SPZLINE'].data['LINEZ'][15])

minw = fits_file['SPECOBJ'].data['WAVEMIN']
maxw = fits_file['SPECOBJ'].data['WAVEMAX']

data_n = len(flux_values)

wave = 10**fits_file['COADD'].data['loglam'] #wavelegnths corresponding to each flux value

#plot
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

fig.savefig('spectrum_0.png')
plt.close('all')