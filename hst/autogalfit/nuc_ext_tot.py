from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os


galaxies =              ['J0826','J0901','J0944','J1107','J1219','J1558','J2116']
nuc_best_mass =  np.array([10.27, 10.11,  10.20,  10.11,  10.45,   9.85,   10.68])
nuc_lo_mass =     np.array([0.04,  0.06,   0.16,   0.32,   0.07,   0.11,    0.26])
nuc_up_mass =     np.array([0.05,  0.05,   0.11,   0.33,   0.10,   0.11,    0.15])

ext_best_mass =   np.array([10.79, 10.69,  10.56, 10.74,  11.03,  10.77,   10.98])
ext_lo_mass =     np.array([0.13,  0.12,   0.12,   0.08,   0.09,   0.07,    0.11])
ext_up_mass =     np.array([0.15,  0.12,   0.13,   0.10,   0.09,   0.07,    0.14])

dew_best_mass =  np.array([10.90, 10.81,  10.80,  10.89,   11.11, 10.77,   11.11])
dew_lo_mass =    np.array([ 0.03,  0.03,   0.05,   0.04,    0.05,  0.05,    0.07])
dew_up_mass =    np.array([ 0.06,  0.05,   0.06,   0.04,    0.06,  0.06,    0.09])


name = 'nuc_ext_tot.pdf'


numbers = np.arange(110)/100.

with PdfPages(name) as pdf:

    fig = plt.figure()

    yo = np.arange(101)/50+9.5

    nucext = np.log10(10**(ext_best_mass) + 10**(nuc_best_mass))

    add_small = np.log10(10**(ext_best_mass-ext_lo_mass) + 10**(nuc_best_mass-nuc_lo_mass))
    add_big = np.log10(10**(ext_best_mass+ext_lo_mass) + 10**(nuc_best_mass+nuc_lo_mass))

    add_err_lo = np.sqrt((nucext - add_small)**2+0.1**2)
    add_err_hi = np.sqrt((add_big - nucext)**2+0.1**2)
    
    tot_err_lo = np.sqrt(dew_lo_mass**2+0.1**2)
    tot_err_hi = np.sqrt(dew_up_mass**2+0.1**2)


    plt.scatter(dew_best_mass, nucext)
    plt.errorbar(dew_best_mass, nucext, yerr=[add_err_lo,add_err_hi], xerr=[tot_err_lo, tot_err_hi], fmt='o', elinewidth=1)

    plt.plot(yo,yo)
    plt.plot(yo,yo, label=r'$\mathcal{M}_{*,central}+\mathcal{M}_{*,extended}=\mathcal{M}_{*,total}$' )
    
    plt.xlabel(r'$\log_{10}(\mathcal{M}_{*, total\ }/\mathcal{M}_{\odot})$', fontsize=14)
    plt.ylabel(r'$\log_{10}((\mathcal{M}_{*, central\ }+\mathcal{M}_{*, extended\ })/\mathcal{M}_{\odot})$', fontsize=14)

    plt.legend(loc='upper left', fontsize=14)

    plt.xlim(10.6,11.4)
    plt.ylim(10.5,11.5)

    #for i in range(0,len(ext_best_mass)):
    #    plt.text(dew_best_mass[i], nucext[i], galaxies[i])
    
    pdf.savefig()
    plt.close()

    #fig = plt.figure()
    #
    #plt.scatter(dew_best_mass, dew_best_mass - nucext)
    #
    #plt.errorbar(dew_best_mass, dew_best_mass - nucext, yerr=[np.sqrt(add_err_lo**2+tot_err_lo**2), np.sqrt(add_err_hi**2+tot_err_hi**2)], xerr=[tot_err_lo, tot_err_hi], fmt='o', elinewidth=1)
    #
    #
    #plt.ylim([-0.3,0.3])
    #plt.xlim([10.6,11.3])
    # 
    #pdf.savefig()
    #plt.close()

# open the pdf file
os.system('open %s &' % name)
