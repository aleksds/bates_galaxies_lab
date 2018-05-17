# Aleks Diamond-Stanic, May 16, 2018
# goal: start with color_mag_plot.py, make code that compares magnitudes, colors, and sizes for two runs of GALFIT
# requires names of two directories as input, current only works for 'independent' runs of galfit
# run color_mag_comparison.py 20180511-1229_psf_independent 20180511-1229_sersic_independent
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

one = sys.argv[1]
two = sys.argv[2]

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']

redshifts = [0.603,0.459,0.712,0.514,0.467,0.451,0.451,0.658,0.608,0.402,0.449,0.728,0.752]
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

model = [one,two]
mags = np.zeros([2,12,2]) #order: psf models then sersic models, normal order for galaxies and filters
chi = np.zeros([2,12,2])
sizepix = np.zeros([12,2])
for m in range(0,2):
    for w in range(0,12):
        for i in range(0,2):
            file = model[m]+'/'+galaxies[w]+'_'+filters[i]+'_fine.galfit.01.band'
            with open(file) as f:
                content = f.readlines()
            mags[m][w][i] = np.float(content[47][4:10])
            chi[m][w][i] = np.float(content[3][14:19])
            sizepix[w][i] = np.float(content[48][4:8])
            
#color magnitude plot
one_mag_475 = np.zeros(12)
one_mag_814 = np.zeros(12)
one_color = np.zeros(12)
two_mag_475 = np.zeros(12)
two_mag_814 = np.zeros(12)
two_color = np.zeros(12)
for w in range(0,12):
    one_mag_475[w] = mags[0][w][0]
    one_mag_814[w] = mags[0][w][1]
    one_color[w] = mags[0][w][0] - mags[0][w][1]
    two_mag_475[w] = mags[1][w][0]
    two_mag_814[w] = mags[1][w][1]    
    two_color[w] = mags[1][w][0] - mags[1][w][1]

name_cm = 'color_vs_mag_'+one+'_'+two+'.pdf'
with PdfPages(name_cm) as pdf:   
    plt.figure()
    plt.scatter(one_mag_814,one_color, label=one, marker='o', color='orange')
    plt.scatter(two_mag_814,two_color, label=two, marker='^', color='purple')
    plt.xlabel('mag_F814W')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Color Magnitude Plot')
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

name_cc = 'color_vs_chi_'+one+'_'+two+'.pdf'
with PdfPages(name_cc) as pdf:   
    plt.figure()
    
    plt.scatter(psfxvals_four,psfyvals, label=one+' (F475W chi)', marker='o', color='blue')
    plt.scatter(sersicxvals_four,sersicyvals, label=two+' (F475W chi)', marker='^', color='blue')
    plt.scatter(psfxvals_eight,psfyvals, label=one+' (F814W chi)', marker='o', color='green')
    plt.scatter(sersicxvals_eight,sersicyvals, label=two+' (F814W chi)', marker='^', color='green')

    plt.xlabel('chi-squared/nu')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Color vs chi-squared/nu')
    plt.legend(loc='upper right')
    pdf.savefig()
    plt.close()
os.system('open %s &' % name_cc)

#color vs size for sersic fits
size_color = np.zeros(12)
size_four = np.zeros(12)
size_eight = np.zeros(12)

kpcrad=np.zeros([12,2])
for w in range(0,12):
    arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshifts[w])
    for i in range(0,2):
        kpcrad[w][i] = (0.025*sizepix[w][i])/arcsecperkpc.value

for w in range(0,12):
    size_four[w] = kpcrad[w][0]
    size_color[w] = mags[1][w][0] - mags[1][w][1]
    size_eight[w] = kpcrad[w][1] 

name_cs = 'color_v_size_'+one+'_'+two+'.pdf'
with PdfPages(name_cs) as pdf:   
    plt.figure()
    
    plt.scatter(size_four,size_color, label='F475W size', marker='^', color='blue')
    plt.scatter(size_eight,size_color, label='F814W size', marker='^', color='green')
    
    plt.xlabel('size (kpc)')
    plt.ylabel('magF475W - magF814W')
    plt.title('Color vs Size')
    plt.legend(loc='lower right')
    pdf.savefig()
    plt.close()
os.system('open %s &' % name_cs)

#comparison plots of mags, colors, chis, and sizes
name_co = 'comparison_'+one+'_'+two+'.pdf'
with PdfPages(name_co) as pdf:   
    fig = plt.figure()

    plt.scatter(one_mag_475, one_mag_475-two_mag_475, marker='o', color='blue')
    plt.xlabel('magF475W_'+one)
    plt.ylabel('magF475W - magF475W')    
    plt.title('magF475W comparison')
    
    pdf.savefig()
    plt.close

    fig = plt.figure()
    
    plt.scatter(one_mag_814, one_mag_814-two_mag_814, marker='o', color='green')
    plt.xlabel('magF814W_'+one)
    plt.ylabel('magF814W - magF814W')    
    plt.title('magF814W comparison')
    
    pdf.savefig()
    plt.close

    fig = plt.figure()
    
    plt.scatter(two_mag_475, one_mag_475-two_mag_475, marker='o', color='red')
    plt.xlabel('magF475W_'+two)
    plt.ylabel('magF475W - magF475W')
    plt.title('magF475W comparison 2')
    
    pdf.savefig()
    plt.close

    fig = plt.figure()
    
    plt.scatter(two_mag_814, one_mag_814-two_mag_814, marker='o', color='purple')
    plt.xlabel('magF814W_'+two)
    plt.ylabel('magF814W - magF814W')
    plt.title('magF814W comparison 2')
    
    pdf.savefig()
    plt.close

    fig = plt.figure()

    plt.scatter(one_color, one_color-two_color, marker='o', color='orange')
    plt.xlabel('color_'+one)
    plt.ylabel('difference in color')
    plt.title('color comparison')

    pdf.savefig()
    plt.close

    fig = plt.figure()

    plt.scatter(psfxvals_four, psfxvals_four-sersicxvals_four, label='F475W', marker='o', color='blue')
    plt.scatter(psfxvals_eight, psfxvals_eight-sersicxvals_eight, label='F814W', marker='o', color='green')
    plt.xlabel('psf chi square/nu values')
    plt.ylabel('difference in chi sqr/nu values between psf and sersic')
    plt.title('chi sqr/nu comparison')
    plt.legend(loc='lower right')

    pdf.savefig()
    plt.close

    fig = plt.figure()

    plt.scatter(size_four, size_four-size_eight, marker='o', color='red')
    plt.xlabel('size_four(kpc)')
    plt.ylabel('size_four-size_eight')
    plt.title('size comparison')

    pdf.savefig()
    plt.close()
    
os.system('open %s &' % name_co)
