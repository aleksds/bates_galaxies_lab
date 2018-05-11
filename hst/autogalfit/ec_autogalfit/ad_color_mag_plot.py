#code to read galfit output files and extract the desired information (integrated magnitudes)
# ad20180511 -- modified to work with files in ec_autogalfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

date_time = sys.argv[1]
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']

redshifts = [0.603,0.459,0.712,0.514,0.467,0.451,0.451,0.658,0.608,0.402,0.449,0.728,0.752]
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

model = ['psf','sersic']
mags = np.zeros([2,12,2]) #order: psf models then sersic models, normal order for galaxies and filters
chi = np.zeros([2,12,2])
sizepix = np.zeros([12,2])
for m in range(0,2):
    for w in range(0,12):
        for i in range(0,2):
            file = date_time+'_'+model[m]+'_independent/'+galaxies[w]+'_'+filters[i]+'_'+model[m]+'_output.galfit.01.band'
            with open(file) as f:
                content = f.readlines()
            mags[m][w][i] = np.float(content[47][4:10])
            chi[m][w][i] = np.float(content[3][14:19])
            sizepix[w][i] = np.float(content[48][4:8])
            
#color magnitude plot
psfyvals = np.zeros(12)
psfxvals = np.zeros(12)
sersicyvals = np.zeros(12)
sersicxvals = np.zeros(12)
for w in range(0,12):
    psfxvals[w] = mags[0][w][1]
    sersicxvals[w] = mags[1][w][1]
    psfyvals[w] = mags[0][w][0] - mags[0][w][1]
    sersicyvals[w] = mags[1][w][0] - mags[1][w][1]


name_cm = 'psf_sersic_color_vs_mag_'+date_time+'.pdf'
with PdfPages(name_cm) as pdf:   
    fig = plt.figure()
    plt.scatter(psfxvals,psfyvals, label='PSF')
    plt.scatter(sersicxvals,sersicyvals, label='SERSIC')
    plt.xlabel('mag_F814W')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Color Magnitude Plot (with psf and sersic fits)')
    plt.legend(loc='upper right')
    pdf.savefig()
    plt.close()
os.system('open %s &' % name_cm)
    
#color vs chi-squared plot
psfyvals = np.zeros(12)
psfxvals_four = np.zeros(12)
sersicyvals = np.zeros(12)
sersicxvals_four = np.zeros(12)
psfxvals_eight = np.zeros(12)
sersicxvals_eight = np.zeros(12)

for w in range(0,12):
    psfxvals_four[w] = chi[0][w][0]
    sersicxvals_four[w] = chi[1][w][0]
    psfyvals[w] = mags[0][w][0] - mags[0][w][1]
    sersicyvals[w] = mags[1][w][0] - mags[1][w][1]
    psfxvals_eight[w] = chi[0][w][1]
    sersicxvals_eight[w] = chi[1][w][1]

name_cc = 'psf_sersic_color_vs_chi_'+date_time+'.pdf'
with PdfPages(name_cc) as pdf:   
    fig = plt.figure()
    
    plt.scatter(psfxvals_four,psfyvals, label='PSF (F475W chi)')
    plt.scatter(sersicxvals_four,sersicyvals, label='SERSIC (F475W chi)')
    plt.scatter(psfxvals_eight,psfyvals, label='PSF (F814W chi)')
    plt.scatter(sersicxvals_eight,sersicyvals, label='SERSIC (F814W chi)')

    plt.xlabel('chi-squared/nu')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Color Magnitude Plot (with psf and sersic fits)')
    plt.legend(loc='upper right')
    pdf.savefig()
    plt.close()
os.system('open %s &' % name_cc)
    
#color vs size for sersic fits
sersicyvals = np.zeros(12)
sersicxvals_four = np.zeros(12)
sersicxvals_eight = np.zeros(12)

kpcrad=np.zeros([12,2])
for w in range(0,12):
    arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshifts[w])
    for i in range(0,2):
        kpcrad[w][i] = (0.025*sizepix[w][i])/arcsecperkpc.value

for w in range(0,12):
    sersicxvals_four[w] = kpcrad[w][0]
    sersicyvals[w] = mags[1][w][0] - mags[1][w][1]
    sersicxvals_eight[w] = kpcrad[w][1] 

name_cs = 'sersic_color_v_size_'+date_time+'.pdf'
with PdfPages(name_cs) as pdf:   
    fig = plt.figure()
    
    plt.scatter(sersicxvals_four,sersicyvals, label='F475W size')
    plt.scatter(sersicxvals_eight,sersicyvals, label='F814W size')
    
    plt.xlabel('size (kpc)')
    plt.ylabel('magF475W - magF814W')
    plt.title('Color vs Size')
    plt.legend(loc='lower right')
    pdf.savefig()
    plt.close()
os.system('open %s &' % name_cs)
