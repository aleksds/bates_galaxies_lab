
import numpy as np
import matplotlib.pyplot as plt
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


log_mass = np.arange(100)/33+9.

log_sig1_quie_260 = (log_mass - 10.5) * 0.67 + 9.80
log_sig1_sf_260 = (log_mass - 10.5) * 0.89 + 9.33

log_sig1_quie_180 = (log_mass - 10.5) * 0.64 + 9.74
log_sig1_sf_180 = (log_mass - 10.5) * 0.86 + 9.25

log_sig1_quie_120 = (log_mass - 10.5) * 0.65 + 9.64
log_sig1_sf_120 = (log_mass - 10.5) * 0.88 + 9.16

log_sig1_quie_075 = (log_mass - 10.5) * 0.65 + 9.53
log_sig1_sf_075 = (log_mass - 10.5) * 0.89 + 9.12

log_sige_quie_260 = (log_mass - 10.5) * (-0.52) + 10.28
log_sige_sf_260 = (log_mass - 10.5) * 0.56 + 8.8

log_sige_quie_180 = (log_mass - 10.5) * (-0.52) + 9.91
log_sige_sf_180 = (log_mass - 10.5) * 0.64 + 8.68

log_sige_quie_120 = (log_mass - 10.5) * (-0.45) + 9.53
log_sige_sf_120 = (log_mass - 10.5) * 0.60 + 8.54

log_sige_quie_075 = (log_mass - 10.5) * (-0.42) + 9.19
log_sige_sf_075 = (log_mass - 10.5) * 0.60 + 8.46


filename = 'sigma_mass.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    ax = fig.add_subplot(2,2,2)

    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')

    ax.plot(10**(log_mass), 10**(log_sig1_quie_260), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_260), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(3e8,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(1.1e10,2e10,'2.2<z<3.0')    

    ax = fig.add_subplot(2,2,1)

    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')

    ax.plot(10**(log_mass), 10**(log_sig1_quie_075), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_075), color='blue', linestyle='dashed')

    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(3e8,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(1.1e10,2e10,'0.5<z<1.0')       
    
    pdf.savefig()
    plt.close()


#os.system('open %s &' % filename)    
    
filename = 'sigma_1_mass_4panels.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    # first plot: 2.2<z<3.0
    ax = fig.add_subplot(2,2,1)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')
    #plt.scatter(10**tot_best_mass, sigma_star_best, marker='o')
    #plt.scatter(10**nuc_best_mass, sigma_star_1kpc, marker='s')

    ax.plot(10**(log_mass), 10**(log_sig1_quie_260), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_260), color='blue', linestyle='dashed')
    
    #log_sige_quie_075 = (log_mass - 10.5) * -0.42 + 9.19
    #plt.plot(10**(log_mass), 10**(log_sige_quie_075), color='red')

    #log_sige_quie_260 = (log_mass - 10.5) * -0.52 + 10.28
    #plt.plot(10**(log_mass), 10**(log_sige_quie_260), color='red')

    #log_sig1_quie_075 = (log_mass - 10.5) * 0.65 + 9.53
    #plt.plot(10**(log_mass), 10**(log_sig1_quie_075), color='red', linestyle='dashed')
    
    ax.plot(10**(log_mass), 10**(log_sig1_quie_260), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_260), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'2.2<z<3.0')


    # second plot: 1.4<z<2.2
    ax = fig.add_subplot(2,2,2)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_180), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_180), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'1.4<z<2.2')

    # third plot: 1.0<z<1.4
    ax = fig.add_subplot(2,2,3)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_120), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_120), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'1.0<z<1.4')

    # fourth plot: 0.5<z<1.0
    ax = fig.add_subplot(2,2,4)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_075), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_075), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'0.5<z<1.0')
    
    #plt.legend(fontsize=12, loc='lower left')

    #for i in range(0, len(galaxies)):
    #    plt.text(sigma_star_best[i], vflow[i], galaxies[i])
    
    pdf.savefig()
    plt.close()

    

#os.system('open %s &' % filename)


filename = 'sigma_e_mass_4panels.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    # first plot: 2.2<z<3.0
    ax = fig.add_subplot(2,2,1)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*')

    ax.plot(10**(log_mass), 10**(log_sige_quie_260), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_260), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,1e12)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'2.2<z<3.0')


    # second plot: 1.4<z<2.2
    ax = fig.add_subplot(2,2,2)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*')
    ax.plot(10**(log_mass), 10**(log_sige_quie_180), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_180), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,1e12)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'1.4<z<2.2')

    # third plot: 1.0<z<1.4
    ax = fig.add_subplot(2,2,3)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*')
    ax.plot(10**(log_mass), 10**(log_sige_quie_120), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_120), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,1e12)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'1.0<z<1.4')

    # fourth plot: 0.5<z<1.0
    ax = fig.add_subplot(2,2,4)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*')
    ax.plot(10**(log_mass), 10**(log_sige_quie_075), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_075), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(3e7,1e12)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.text(5e10,1e8,'0.5<z<1.0')
    
    #plt.legend(fontsize=12, loc='lower left')

    #for i in range(0, len(galaxies)):
    #    plt.text(sigma_star_best[i], vflow[i], galaxies[i])
    
    pdf.savefig()
    plt.close()

    

#os.system('open %s &' % filename)


filename = 'sigmae_sigma1.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    # first plot: 2.2<z<3.0
    ax = fig.add_subplot(2,2,1)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*', color='#2ca02c')

    ax.plot(10**(log_mass), 10**(log_sige_quie_260), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_260), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(1e8,1e12)
    #ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_e$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    #plt.text(5e10,1e8,'2.2<z<3.0')
    plt.text(1.05e10, 1.6e10, '2.2<z<3.0 quiescent', rotation=-9, fontsize=9)
    plt.text(1.05e10, 6.8e8, '2.2<z<3.0 star-forming', rotation=10, fontsize=9)

    # second plot: 0.5<z<1.0
    ax = fig.add_subplot(2,2,2)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*', color='#2ca02c')
    ax.plot(10**(log_mass), 10**(log_sige_quie_075), color='red', linestyle='solid')
    ax.plot(10**(log_mass), 10**(log_sige_sf_075), color='blue', linestyle='dashed')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(1e8,1e12)
    #ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    #plt.text(1.1e10,5e11,'0.5<z<1.0', fontsize=9)
    plt.text(1.05e10, 1.2e9, '0.5<z<1.0', rotation=-7, fontsize=9)

    # third plot: 2.2<z<3.0
    ax = fig.add_subplot(2,2,3)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*', color='#2ca02c', label='Compact starbursts (this paper)')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_260), color='red', linestyle='solid', label='Quiescent galaxies (Barro+17)')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_260), color='blue', linestyle='dashed', label='Star-forming galaxies (Barro+17)')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(1e8,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    #plt.text(1.1e10,1.8e10,'2.2<z<3.0', fontsize=9)
    plt.text(1.05e10, 2.5e9, '2.2<z<3.0', rotation=24, fontsize=9)

    plt.legend(fontsize=7.5)

    # fourth plot: 0.5<z<1.0
    ax = fig.add_subplot(2,2,4)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*', color='#2ca02c', label='Compact starbursts (z~0.6)')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_075), color='red', linestyle='solid', label='Quiescent galaxies (Barro+17)')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_075), color='blue', linestyle='dashed', label='Star-forming galaxies (Barro+17)')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(1e8,3e10)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    #ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    #plt.text(1.1e10,1.8e10,'0.5<z<1.0', fontsize=9)
    plt.text(1.05e10, 1.6e9, '0.5<z<1.0', rotation=24, fontsize=9)


    #plt.legend(fontsize=7)
    
    pdf.savefig()
    plt.close()

    

os.system('open %s &' % filename)


filename = 'sigmae_sigma1_combined.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    # first plot: Sigma_e
    ax = fig.add_subplot(2,1,1)
    
    ax.scatter(10**tot_best_mass, sigma_star_best, marker='*', color='#2ca02c', label='Compact starbursts 0.4<z<0.8')

    ax.plot(10**(log_mass), 10**(log_sige_quie_260), color='red', linestyle='solid', label='Quiescent galaxies 2.2<z<3.0')
    ax.plot(10**(log_mass), 10**(log_sige_sf_260), color='blue', linestyle='dashed', label='Star-forming galaxies 2.2<z<3.0')

    ax.plot(10**(log_mass), 10**(log_sige_quie_075), color='red', linestyle='dashdot', label='Quiescent galaxies 0.5<z<1.0')
    ax.plot(10**(log_mass), 10**(log_sige_sf_075), color='blue', linestyle='dotted', label='Star-forming galaxies 0.5<z<1.0')
    
    ax.set_xlim(1e10, 5e11)
    ax.set_ylim(1e8,1e12)
    #ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_e$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.text(5e10,1e8,'2.2<z<3.0')
    #plt.legend()

    # second plot: Sigma_1
    ax = fig.add_subplot(2,1,2)
    
    ax.scatter(10**tot_best_mass, sigma_star_1kpc, marker='*', color='#2ca02c', label='Compact starbursts 0.4<z<0.8')
    ax.plot(10**(log_mass), 10**(log_sig1_quie_260), color='red', linestyle='solid', label='Quiescent galaxies 2.2<z<3.0')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_260), color='blue', linestyle='dashed', label='Star-forming galaxies 2.2<z<3.0')

    ax.plot(10**(log_mass), 10**(log_sig1_quie_075), color='red', linestyle='dashdot', label='Quiescent galaxies 0.5<z<1.0')
    ax.plot(10**(log_mass), 10**(log_sig1_sf_075), color='blue', linestyle='dotted', label='Star-forming galaxies 0.5<z<1.0')
    
    ax.set_xlim(1e9, 5e11)
    ax.set_ylim(1e8,1e12)
    ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ax.set_ylabel(r'$\Sigma_1$ [M$_\odot$ kpc$^{-2}$]', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.text(5e10,1e8,'2.2<z<3.0')

    plt.legend(fontsize=8)
    
    pdf.savefig()
    plt.close()

    

#os.system('open %s &' % filename)
