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
nuc_best_mass =  np.array([10.27, 10.11,  10.41,  10.20,  10.11,  10.45,  10.51,  10.18,  9.85,  10.65,  10.68,  10.56])
nuc_lo_mass =     np.array([0.04,  0.06,   0.24,   0.16,   0.32,   0.07,   0.25,   0.20,  0.11,   0.13,   0.26,   0.07])
nuc_up_mass =     np.array([0.05,  0.05,   0.29,   0.11,   0.33,   0.10,   0.15,   0.21,  0.11,   0.07,   0.15,   0.10])

nuc_lo_mass = np.sqrt(nuc_lo_mass**2 + 0.1**2)
nuc_up_mass = np.sqrt(nuc_up_mass**2 + 0.1**2)

nuc_mass_loval = 10**(nuc_best_mass - nuc_lo_mass)
nuc_mass_hival = 10**(nuc_best_mass + nuc_up_mass)

tot_best_mass =  np.array([10.90, 10.81,  10.98,  10.80,  10.89,   11.11, 10.86,  10.84,  10.77,  11.13,  11.11,  11.16])
tot_lo_mass =    np.array([ 0.03,  0.03,   0.03,   0.05,   0.04,    0.05,  0.02,   0.04,   0.05,   0.04,   0.07,   0.06])
tot_up_mass =    np.array([ 0.06,  0.05,   0.05,   0.06,   0.04,    0.06,  0.04,   0.05,   0.06,   0.05,   0.09,   0.05])

tot_lo_mass = np.sqrt(tot_lo_mass**2 + 0.1**2)
tot_up_mass = np.sqrt(tot_up_mass**2 + 0.1**2)

tot_mass_loval = 10**(tot_best_mass - tot_lo_mass)
tot_mass_hival = 10**(tot_best_mass + tot_up_mass)

re_best = np.array([0.0151, 0.0149, 0.0105, 0.0099, 0.0156, 0.0257, 0.0127, 0.0118, 0.0387, 0.1289, 0.0216, 0.0145]) 
re_unc = np.array([0.0031, 0.0033, 0.0027, 0.0030, 0.0041, 0.0038, 0.0023, 0.0025, 0.0064, 0.0157, 0.0046, 0.0045])
re_best_arc = re_best * u.arcsec

z = np.array([0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752])


cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

re_best_kpc = re_best_arc / cosmo.arcsec_per_kpc_proper(z)
re_best_kpc_hi = (re_best_arc + re_unc * u.arcsec) / cosmo.arcsec_per_kpc_proper(z)
re_best_kpc_lo = (re_best_arc - re_unc * u.arcsec) / cosmo.arcsec_per_kpc_proper(z)

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
log_mass_van = log_mass_quie + 0.29897
log_re_van = log_mass_van - 10.7

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

    ax.scatter(10**tot_best_mass, re_best_kpc, marker='.', color='orange')
    ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*', color='green')
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

    ax.scatter(10**tot_best_mass, re_best_kpc, marker='.', color='orange')
    ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*', color='green')
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

#os.system('open %s &' % filename)    
    
filename = 'size_mass_onepanel.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()


    ax = fig.add_subplot(1,1,1)

    ax.scatter(10**tot_best_mass, re_best_kpc, marker='.', color='#ff7f0e', label=r'$\mathcal{M}_{*,total}$')#label='compact starbursts, M*=Mtot')
    ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*', color='#ff7f0e', label=r'$\mathcal{M}_{*,central}$')#label='compact starbursts, M*=Mnuc') # label=r'compact starbursts: $\mathcal{M}_{*,central}$')
    ax.plot(10**log_mass_quie, re_early_075, color='red', label='early-type galaxies', linestyle='dashdot')#linestyle=(0, (5, 3))) # densely dashed
    ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle=(0, (5, 5))) # loosely dashed
    ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle=(0, (5, 5))) # loosely dashed

    #ax.errorbar(10**nuc_best_mass, np.array(re_best_kpc), yerr=[np.array(re_best_kpc/u.kpc) - np.array(re_best_kpc_lo/u.kpc), np.array(re_best_kpc_hi/u.kpc) - np.array(re_best_kpc/u.kpc)], xerr=[10**nuc_best_mass-nuc_mass_loval,nuc_mass_hival-10**nuc_best_mass], fmt='none', color='green', elinewidth=1)
    
    ax.errorbar(10**tot_best_mass, np.array(re_best_kpc), yerr=[np.array(re_best_kpc/u.kpc) - np.array(re_best_kpc_lo/u.kpc), np.array(re_best_kpc_hi/u.kpc) - np.array(re_best_kpc/u.kpc)], xerr=[10**tot_best_mass-tot_mass_loval,tot_mass_hival-10**tot_best_mass], fmt='none', color='#ff7f0e', elinewidth=1)

    eb = ax.errorbar(10**tot_best_mass, np.array(re_best_kpc), xerr=[10**tot_best_mass-10**nuc_best_mass, np.zeros(len(tot_best_mass))], fmt='none', elinewidth=1, color='#ff7f0e')
    eb[-1][0].set_linestyle((0, (1, 5))) #eb1[-1][0] is the LineCollection objects of the errorbar lines
    
    #ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dotted')
    #ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dashed')

    
    #ax.plot(10**log_mass_quie, re_early_075, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle='dashed')

    ax.set_xlim(4e9, 4e11)
    ax.set_ylim(0.04,20)
    ax.set_xlabel(r'$\mathcal{M}_*$ [$\mathcal{M}_\odot$]', fontsize=13)
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
    
    #ax = fig.add_subplot(2,2,2)

    #ax.scatter(10**tot_best_mass, re_best_kpc, marker='.', color='orange')
    #ax.scatter(10**nuc_best_mass, re_best_kpc, marker='*', color='green')
    ax.plot(10**log_mass_quie, re_early_275, color='red', linestyle='dashdot')#(0, (5, 1))) # densely dashed
    ax.plot(10**log_mass_quie, re_early_275_lo, color='red', linestyle=(0, (5, 5))) # loosely dashed
    ax.plot(10**log_mass_quie, re_early_275_hi, color='red', linestyle=(0, (5, 5))) # loosely dashed


    ax.plot(10**log_mass_van, 10**log_re_van, color='black', linestyle='solid', label='compact massive galaxies')

    log_re_tmp = np.arange(20)/10.-2
    log_mass_tmp = np.zeros(len(log_re_tmp))+10.6
    ax.plot(10**log_mass_tmp, 10.**log_re_tmp, color='black', linestyle='solid')#linestyle='dashdot')
    
    
    plt.text(9.4e9, 0.37, '2.5<z<3.0', rotation=23, fontsize=11)
    plt.text(9.4e9, 1.8, '0.5<z<1.0', rotation=22, fontsize=11)

    
    plt.legend(loc='upper left')
    
    #ax.plot(10**log_mass, re_late_275, color='blue', linestyle='dotted')
    #ax.plot(10**log_mass, re_late_075, color='blue', linestyle='dotted')

    
    #ax.plot(10**log_mass_quie, re_early_075, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_lo, color='red', linestyle='dashed')
    #ax.plot(10**log_mass_quie, re_early_075_hi, color='red', linestyle='dashed')

    #ax.set_xlim(1e9, 5e11)
    #ax.set_ylim(0.05,20)
    #ax.set_xlabel(r'$M\star$ [M$_\odot$]', fontsize=13)
    ##ax.set_ylabel(r'$r_e$ [kpc]', fontsize=13)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.tick_params(axis='both', which='major', labelsize=12)
    #plt.text(1.1e10,2e10,'2.2<z<3.0')    

    #ax.set_yticks([0.1, 1, 10])
    #ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    #labels = [item.get_text() for item in ax.get_yticklabels()]
    #print(labels)
    #labels[0] = '0.1'
    #labels[1] = '1'
    #labels[2] = '10'

    #ax.set_yticklabels(labels)
    
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)    
