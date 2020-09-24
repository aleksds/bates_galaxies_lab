
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
#vflow = np.array([           1600,  1700,   3000,   2100,   1828,   2250,   2000,   2050,  1350,   2600,   1900,   1100])
#vmax = np.array([         1456.2, 1575.4, 2884.8, 1878.3, 2015.6, 1962.6, 1995.2, 1768.8,1151.2, 2374.1, 1915.1,  950.4])
#vavg = np.array([         1131.2, 1210.1, 2395.3, 1172.4, 1359.8, 1543.7,  766.9, 1378.9, 807.1, 1437.5, 1056.2,  342.8])

vavg = np.array(   [1247.52,1308.73,2514.93,1286.62,1416.76,1684.03, 760.78,1288.30,868.45, 807.53,1124.12, 512.42])
vmax = np.array(   [1657.28,1777.88,3019.30,1988.41,2050.23,2031.88,1997.66,2447.93,1255.89,2503.18,2238.64,1042.26])


nuc_up_mass = np.sqrt(nuc_up_mass**2 + 0.1**2)
nuc_lo_mass = np.sqrt(nuc_lo_mass**2 + 0.1**2)

tot_best_mass = np.array([10.57, 10.47, 10.67, 10.52, 10.57, 11.40, 10.53, 10.57, 10.79, 11.50, 10.94, 10.92])
tot_up_mass = np.array([0.01, 0.01, 0.03, 0.03, 0.02, 0.11, 0.04, 0.02, 0.09, 0.10, 0.06, 0.02])
tot_up_mass = np.sqrt(tot_up_mass**2 + 0.1**2)
tot_lo_mass = np.array([0.01, 0.01, 0.03, 0.02, 0.01, 0.11, 0.03, 0.02, 0.06, 0.11, 0.06, 0.03])
tot_lo_mass = np.sqrt(tot_lo_mass**2 + 0.1**2)

re_best = np.array([0.0151, 0.0149, 0.0105, 0.0099, 0.0156, 0.0257, 0.0127, 0.0118, 0.0387, 0.1289, 0.0216, 0.0145]) 
re_unc = np.array([0.0031, 0.0033, 0.0027, 0.0030, 0.0041, 0.0038, 0.0023, 0.0025, 0.0064, 0.0157, 0.0046, 0.0045])
re_best_arc = re_best * u.arcsec

z = np.array([0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752])


cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

re_best_kpc = re_best_arc / cosmo.arcsec_per_kpc_proper(z)
vesc_best = np.sqrt(G * 10**nuc_best_mass * const.M_sun / re_best_kpc).to('km/s') * u.s / u.km
vesc_hi = np.sqrt(G * 10**(nuc_best_mass + nuc_up_mass) * const.M_sun / (re_best_kpc-(re_unc/re_best*re_best_kpc))).to('km/s') * u.s / u.km
vesc_lo = np.sqrt(G * 10**(nuc_best_mass - nuc_lo_mass) * const.M_sun / (re_best_kpc+(re_unc/re_best*re_best_kpc))).to('km/s') * u.s / u.km

vesc_hi_re = np.sqrt(G * 10**(nuc_best_mass) * const.M_sun / (re_best_kpc-(re_unc/re_best*re_best_kpc))).to('km/s') * u.s / u.km
vesc_lo_re = np.sqrt(G * 10**(nuc_best_mass) * const.M_sun / (re_best_kpc+(re_unc/re_best*re_best_kpc))).to('km/s') * u.s / u.km

vesc_hi_mass = np.sqrt(G * 10**(nuc_best_mass + nuc_up_mass) * const.M_sun / re_best_kpc).to('km/s') * u.s / u.km
vesc_lo_mass = np.sqrt(G * 10**(nuc_best_mass - nuc_lo_mass) * const.M_sun / re_best_kpc).to('km/s') * u.s / u.km


print('mean vmax / vesc,central', np.mean(vmax/vesc_best))
print('median vmax / vesc,central', np.median(vmax/vesc_best))
print('mean vavg / vesc,central', np.mean(vavg/vesc_best))
print('median vavg / vesc,central', np.median(vavg/vesc_best))

print('mean vmax', np.mean(vmax))
print('mean vavg', np.mean(vavg))

filename = 'vesc_flip.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    vel_array = np.arange(101)*40

    plt.scatter(vesc_best, vmax, marker='o', color='#ff7f0e', label=r'$v_{max}$')
    eb = plt.errorbar(vesc_best, vmax, yerr=[vmax-vavg,np.zeros(len(vavg))], fmt='none', elinewidth=1, color='#ff7f0e', label=r'[$v_{avg}$, $v_{max}$]')
    eb[-1][0].set_linestyle((0, (1, 3)))#'dotted') #eb1[-1][0] is the LineCollection objects of the errorbar lines
    plt.scatter(vesc_best, vavg, marker='+', color='#ff7f0e', label=r'$v_{avg}$')


    

    plt.errorbar(vesc_best, vmax, xerr=[vesc_best-vesc_lo, vesc_hi-vesc_best], fmt='none', elinewidth=1, label=r'$r_e$ and mass uncertainty', color='#1f77b4')
    #plt.plot(vel_array, vel_array/3, linestyle='--')
    #plt.errorbar(vmax, vesc_best, yerr=[vesc_best-vesc_lo_mass, vesc_hi_mass-vesc_best], fmt='o')

    plt.errorbar(vesc_best, vmax, xerr=[vesc_best-vesc_lo_re, vesc_hi_re-vesc_best], fmt='none', elinewidth=3, label=r'$r_e$ uncertainty', color='#ff7f0e') 

    plt.ylim(0,3200)
    plt.xlim(0,3200)
    plt.ylabel(r'Observed Outflow Velocity [km s$^{-1}$]', fontsize=13)
    plt.xlabel(r'Central Escape Velocity [km s$^{-1}$]', fontsize=13)
    plt.plot(vel_array, vel_array, label=r'$v_{outflow} = v_{escape}$', linestyle='dashed', color='#2ca02c')
    #plt.plot(vel_array, vel_array/3, label=r'$3\times v_{esc} = v_{out}$')

    plt.tick_params(axis='both', which='major', labelsize=12)


#['#1f77b4', # blue
# '#ff7f0e', # orange
# '#2ca02c', # green
# '#d62728', # red
# '#9467bd',
# '#8c564b',
# '#e377c2',
# '#7f7f7f',
# '#bcbd22',
# '#17becf']
    
    plt.legend(fontsize=11, loc='lower right')

    #plt.scatter(vesc_best, vavg)

    #for i in range(0, len(galaxies)):
    #    plt.text(vmax[i], vesc_best[i], galaxies[i])

    plt.tight_layout()
    
    pdf.savefig()
    plt.close()

    
os.system('open %s &' % filename)

filename = 'vesc_avg.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    vel_array = np.arange(101)*40

    plt.scatter(vesc_best, vavg)
    plt.errorbar(vesc_best, vavg, xerr=[vesc_best-vesc_lo, vesc_hi-vesc_best], fmt='.', elinewidth=1, label=r'$r_e$ and mass uncertainty')
    #plt.plot(vel_array, vel_array/3, linestyle='--')
    #plt.errorbar(vavg, vesc_best, yerr=[vesc_best-vesc_lo_mass, vesc_hi_mass-vesc_best], fmt='o')    
    plt.errorbar(vesc_best, vavg, xerr=[vesc_best-vesc_lo_re, vesc_hi_re-vesc_best], fmt='o', elinewidth=3, label=r'$r_e$ uncertainty')
    plt.ylim(0,3200)
    plt.xlim(0,2400)
    plt.ylabel(r'Observed Outflow Velocity [km s$^{-1}$]', fontsize=13)
    plt.xlabel(r'Central Escape Velocity [km s$^{-1}$]', fontsize=13)
    plt.plot(vel_array, vel_array, label=r'$v_{escape} = v_{outflow}$')
    #plt.plot(vel_array, vel_array/3, label=r'$3\times v_{esc} = v_{out}$')
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(fontsize=12, loc='lower right')

    #for i in range(0, len(galaxies)):
    #    plt.text(vavg[i], vesc_best[i], galaxies[i])
    
    pdf.savefig()
    plt.close()

    
#os.system('open %s &' % filename)
