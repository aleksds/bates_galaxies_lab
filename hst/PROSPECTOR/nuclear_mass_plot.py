import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from astropy.io import ascii

penmass = ascii.read('Penultimate_Mass_PPD.csv')

nuc_best_mass = np.array([10.20, 10.07, 10.41, 10.14, 10.11, 10.72, 10.51, 10.18, 9.81, 10.52, 10.68, 10.74])
nuc_lo_mass = np.array([0.24, 0.25, 0.24, 0.18, 0.31, 0.42, 0.25, 0.20, 0.16, 0.23, 0.25, 0.27])
nuc_up_mass = np.array([0.26, 0.28, 0.29, 0.18, 0.34, 0.25, 0.15, 0.21, 0.12, 0.20, 0.15, 0.25])

tot_best_mass = np.array([10.57, 10.47, 10.67, 10.52, 10.57, 11.40, 10.53, 10.57, 10.79, 11.50, 10.94, 10.92])
tot_lo_mass = np.array([0.01, 0.01, 0.03, 0.02, 0.01, 0.11, 0.03, 0.02, 0.06, 0.11, 0.06, 0.03])
tot_lo_mass = np.sqrt(tot_lo_mass**2 + 0.075**2)
tot_up_mass = np.array([0.01, 0.01, 0.03, 0.03, 0.02, 0.11, 0.04, 0.02, 0.09, 0.10, 0.06, 0.02])
tot_up_mass = np.sqrt(tot_up_mass**2 + 0.075**2)

filename = 'nuclear_mass_plot.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    mass_array = np.arange(101)/25.+8.
    
    plt.scatter(tot_best_mass, nuc_best_mass)
    plt.errorbar(x=tot_best_mass, y=nuc_best_mass, xerr=[tot_lo_mass, tot_up_mass],
                 yerr=[nuc_lo_mass,nuc_up_mass], ls='none')
    plt.plot(mass_array, mass_array, linestyle='-', label='Nuclear Stellar Mass = Total Stellar Mass', color='black')
    plt.plot(mass_array, mass_array-(1-np.log10(3)), linestyle=':', label='Nuclear = 30% Total')
    plt.plot(mass_array, mass_array-1, linestyle='--', label='Nuclear = 10% Total')
    plt.ylim(9.6, 11.6)
    plt.xlim(9.6, 11.6)
    plt.xlabel(u'$log_{10}$(Total Stellar Mass [$M_\u2609$])', fontsize=12)
    plt.ylabel(u'$log_{10}$(Nuclear Stellar Mass [$M_\u2609$])', fontsize=12)
    plt.legend()

    for i in range(0, len(penmass)):
        plt.text(tot_best_mass[i], nuc_best_mass[i], penmass['Galaxy'][i])
    
    pdf.savefig()
    plt.close()

    

os.system('open %s &' % filename)
