import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import os
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
from astropy import constants as const
from astropy import units as u

galaxies =              ['J0826','J0901','J0905','J0944','J1107','J1219','J1341','J1506','J1558','J1613','J2116','J2140']
#nuc_best_mass = np.array([10.20, 10.07,  10.41,  10.14,  10.11,  10.72,  10.51,  10.18,  9.81,  10.52,  10.68,  10.74])
nuc_best_mass =  np.array([10.27, 10.11,  10.41,  10.20,  10.11,  10.45,  10.51,  10.18,  9.85,  10.65,  10.68,  10.56])
#nuc_up_mass =    np.array([0.26,  0.28,   0.29,   0.18,   0.34,   0.25,   0.15,   0.21,  0.12,   0.20,   0.15,   0.25])
nuc_up_mass =     np.array([0.04,  0.05,   0.29,   0.11,   0.34,   0.10,   0.15,   0.21,  0.11,   0.07,   0.15,   0.10])
#nuc_lo_mass =    np.array([0.24,  0.25,   0.24,   0.18,   0.31,   0.42,   0.25,   0.20,  0.16,   0.23,   0.25,   0.27])
nuc_lo_mass =     np.array([0.05,  0.06,   0.24,   0.16,   0.31,   0.07,   0.25,   0.20,  0.11,   0.13,   0.25,   0.07])

#vflow = np.array([          1228,  1206,   2470,   1778,   1828,   1830,    875,   1211,   829,   2416,   1456,    606])
vflow = np.array([           1600,  1700,   3000,   2100,   1828,   2250,   2000,   2050,  1350,   2600,   1900,   1100])

nuc_up_mass = np.sqrt(nuc_up_mass**2 + 0.1**2)
nuc_lo_mass = np.sqrt(nuc_lo_mass**2 + 0.1**2)

tot_best_mass = np.array([10.57, 10.47, 10.67, 10.52, 10.57, 11.40, 10.53, 10.57, 10.79, 11.50, 10.94, 10.92])
tot_up_mass = np.array([0.01, 0.01, 0.03, 0.03, 0.02, 0.11, 0.04, 0.02, 0.09, 0.10, 0.06, 0.02])
tot_up_mass = np.sqrt(tot_up_mass**2 + 0.075**2)
tot_lo_mass = np.array([0.01, 0.01, 0.03, 0.02, 0.01, 0.11, 0.03, 0.02, 0.06, 0.11, 0.06, 0.03])
tot_lo_mass = np.sqrt(tot_lo_mass**2 + 0.075**2)

re_best = np.array([0.0151, 0.0149, 0.0105, 0.0099, 0.0156, 0.0257, 0.0127, 0.0118, 0.0387, 0.1289, 0.0216, 0.0145]) 
re_unc = np.array([0.0031, 0.0033, 0.0027, 0.0030, 0.0041, 0.0038, 0.0023, 0.0025, 0.0064, 0.0157, 0.0046, 0.0045])
re_best_arc = re_best * u.arcsec

z = np.array([0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752])


cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

re_best_kpc = re_best_arc / cosmo.arcsec_per_kpc_proper(z)

sigma_star_corr = np.array([0.96/0.5, 0.97/0.5, 0.98/0.5, 0.99/0.5, 0.97/0.5, 0.93/0.5, 0.97/0.5, 0.98/0.5, 0.88/0.5, 0.58/0.5, 0.92/0.5, 0.96/0.5])

sigma_star_1kpc = 10**nuc_best_mass / (2. * np.pi * 1**2) * sigma_star_corr

sigma_star_best = 10**nuc_best_mass / (2. * np.pi * re_best_kpc**2) * u.kpc * u.kpc

sigma_star_hi = 10**(nuc_best_mass + nuc_up_mass) / (2. * np.pi * (re_best_kpc-(re_unc/re_best*re_best_kpc))**2) * u.kpc * u.kpc

sigma_star_lo = 10**(nuc_best_mass - nuc_lo_mass) / (2. * np.pi * (re_best_kpc+(re_unc/re_best*re_best_kpc))**2) * u.kpc * u.kpc

sigma_star_hi_re = 10**(nuc_best_mass) / (2. * np.pi * (re_best_kpc-(re_unc/re_best*re_best_kpc))**2) * u.kpc * u.kpc

sigma_star_lo_re = 10**(nuc_best_mass) / (2. * np.pi * (re_best_kpc+(re_unc/re_best*re_best_kpc))**2) * u.kpc * u.kpc

sigma_star_hi_mass = 10**(nuc_best_mass + nuc_up_mass) / (2. * np.pi * (re_best_kpc)**2) * u.kpc * u.kpc

sigma_star_lo_mass = 10**(nuc_best_mass - nuc_lo_mass) / (2. * np.pi * (re_best_kpc)**2) * u.kpc * u.kpc


#log_mass = np.arange(100)/33+np.log10(3e9)
log_mass = np.arange(67)/33+np.log10(3e9)
#log_mass_quie = np.arange(100)/33+np.log10(2e10)
blah = 100. / (np.log10(3e11) - np.log10(2e10))
log_mass_quie = np.arange(101)/blah+np.log10(2e10)

log_sig1_quie_260 = (log_mass_quie - 10.5) * 0.67 + 9.80
log_sig1_quie_260_lo = (log_mass_quie - 10.5) * 0.67 + 9.80 - 0.2
log_sig1_quie_260_hi = (log_mass_quie - 10.5) * 0.67 + 9.80 + 0.2

log_sig1_sf_260 = (log_mass - 10.5) * 0.89 + 9.33

log_sig1_quie_180 = (log_mass_quie - 10.5) * 0.64 + 9.74
log_sig1_sf_180 = (log_mass - 10.5) * 0.86 + 9.25

log_sig1_quie_120 = (log_mass_quie - 10.5) * 0.65 + 9.64
log_sig1_sf_120 = (log_mass - 10.5) * 0.88 + 9.16

log_sig1_quie_075 = (log_mass_quie - 10.5) * 0.65 + 9.53
log_sig1_quie_075_lo = (log_mass_quie - 10.5) * 0.65 + 9.53 - 0.2
log_sig1_quie_075_hi = (log_mass_quie - 10.5) * 0.65 + 9.53 + 0.2
    
log_sig1_sf_075 = (log_mass - 10.5) * 0.89 + 9.12

log_sige_quie_260 = (log_mass_quie - 10.5) * (-0.52) + 10.28
log_sige_quie_260_lo = (log_mass_quie - 10.5) * (-0.52) + 10.28 - 0.5
log_sige_quie_260_hi = (log_mass_quie - 10.5) * (-0.52) + 10.28 + 0.5
    
log_sige_sf_260 = (log_mass - 10.5) * 0.56 + 8.8

log_sige_quie_180 = (log_mass_quie - 10.5) * (-0.52) + 9.91
log_sige_sf_180 = (log_mass - 10.5) * 0.64 + 8.68

log_sige_quie_120 = (log_mass_quie - 10.5) * (-0.45) + 9.53
log_sige_sf_120 = (log_mass - 10.5) * 0.60 + 8.54

log_sige_quie_075 = (log_mass_quie - 10.5) * (-0.42) + 9.19
log_sige_quie_075_lo = (log_mass_quie - 10.5) * (-0.42) + 9.19 - 0.5
log_sige_quie_075_hi = (log_mass_quie - 10.5) * (-0.42) + 9.19 + 0.5

log_sige_sf_075 = (log_mass - 10.5) * 0.60 + 8.46

re_early_275 = 10.**(-0.06) * (10**(log_mass_quie) / 5e10)**0.79
re_early_275_hi = 10.**(np.log10(re_early_275)+0.19)
re_early_275_lo = 10.**(np.log10(re_early_275)-0.19)

re_late_275 = 10.**(0.51) * (10**(log_mass) / 5e10)**0.18

re_late_075 = 10.**(0.86) * (10**(log_mass) / 5e10)**0.16
    
re_early_075 = 10.**(0.60) * (10**(log_mass_quie) / 5e10)**0.75
re_early_075_hi = 10.**(np.log10(re_early_075)+0.10)
re_early_075_lo = 10.**(np.log10(re_early_075)-0.10)

filename = 'size_mass.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()


    ax = fig.add_subplot(2,2,1)

    ax.scatter(10**tot_best_mass, re_best_kpc, marker='o', color='green')
    ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*')
    ax.plot(10**log_mass_quie, re_early_075, color='red')
    ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle='dashed')
    ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle='dashed')

    ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dotted')
    #ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dashed')

    
    #ax.plot(10**log_mass_quie, re_early_075, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle='dashed')

    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(0.05,20)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$r_e$ [kpc]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.set_yticks([0.1, 1, 10])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    labels = [item.get_text() for item in ax.get_yticklabels()]
    print(labels)
    labels[0] = '0.1'
    labels[1] = '1'
    labels[2] = '10'

    ax.set_yticklabels(labels)
    
    ax = fig.add_subplot(2,2,2)

    ax.scatter(10**tot_best_mass, re_best_kpc, marker='o', color='green')
    ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*')
    ax.plot(10**log_mass_quie, re_early_275, color='red')
    ax.plot(10**log_mass_quie, re_early_275_lo, color='red', linestyle='dashed')
    ax.plot(10**log_mass_quie, re_early_275_hi, color='red', linestyle='dashed')

    ax.plot(10**log_mass, re_late_275, color='blue', linestyle='dotted')
    #ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dashed')

    
    #ax.plot(10**log_mass_quie, re_early_075, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle='dashed')

    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(0.05,20)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$r_e$ [kpc]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    #plt.text(1.1e10,2e10,'2.2<z<3.0')    

    ax.set_yticks([0.1, 1, 10])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    labels = [item.get_text() for item in ax.get_yticklabels()]
    print(labels)
    labels[0] = '0.1'
    labels[1] = '1'
    labels[2] = '10'

    ax.set_yticklabels(labels)
    
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)    
    

