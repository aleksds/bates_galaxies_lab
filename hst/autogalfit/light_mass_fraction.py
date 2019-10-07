from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# function to return a flux in maggies given an AB magnitude
def flux(mag):
    flux = 10. ** (mag / (-2.5))
    return flux

data = ascii.read('bgl_phot.dat')
dtot = ascii.read('bgl_tot.dat')

galaxies =              ['J0826','J0901','J0905','J0944','J1107','J1219','J1341','J1506','J1558','J1613','J2116','J2140']
nuc_best_mass =  np.array([10.27, 10.11,  10.41,  10.20,  10.11,  10.45,  10.51,  10.18,  9.85,  10.65,  10.68,  10.56])
nuc_lo_mass =     np.array([0.04,  0.06,   0.24,   0.16,   0.32,   0.07,   0.25,   0.20,  0.11,   0.13,   0.26,   0.07])
nuc_up_mass =     np.array([0.05,  0.05,   0.29,   0.11,   0.33,   0.10,   0.15,   0.21,  0.11,   0.07,   0.15,   0.10])

ext_best_mass =  np.array([10.66, 10.55,  10.12,  10.54,  10.67,  10.99,  10.43,  11.03,  10.63,   9.88,  10.98,  10.42])
ext_lo_mass =     np.array([0.32,  0.26,   0.18,   0.14,   0.17,   0.16,   0.44,   0.16,   0.15,   0.47,   0.27,   0.25])
ext_up_mass =     np.array([0.21,  0.19,   0.17,   0.13,   0.13,   0.10,   0.43,   0.15,   0.13,   0.47,   0.19,   0.20])

mass_frac = 10**(nuc_best_mass) / (10**(nuc_best_mass) + 10**(ext_best_mass))
mass_frac_hi = 10**(nuc_best_mass+nuc_up_mass) / (10**(nuc_best_mass+nuc_up_mass) + 10**(ext_best_mass-ext_lo_mass))
mass_frac_lo = 10**(nuc_best_mass-nuc_lo_mass) / (10**(nuc_best_mass-nuc_lo_mass) + 10**(ext_best_mass+ext_up_mass))

diff_frac_hi = mass_frac_hi - mass_frac
diff_frac_lo = mass_frac - mass_frac_lo

# F475W
nuc_flux_475 = flux(data['m475'])
nuc_unc_475 = nuc_flux_475 / 1.086 * data['u475']

tot_flux_475 = flux(dtot['m475'])
tot_unc_475 = tot_flux_475 / 1.086 * np.sqrt((dtot['u475'])**2 + 0.03**2)

frac_475 = nuc_flux_475 / tot_flux_475
func_475 = frac_475 * np.sqrt((nuc_unc_475/nuc_flux_475)**2+(tot_unc_475/tot_flux_475)**2)

# F814W
nuc_flux_814 = flux(data['m814'])
nuc_unc_814 = nuc_flux_814 / 1.086 * data['u814']

tot_flux_814 = flux(dtot['m814'])
tot_unc_814 = tot_flux_814 / 1.086 * np.sqrt((dtot['u814'])**2 + 0.03**2)

frac_814 = nuc_flux_814 / tot_flux_814
func_814 = frac_814 * np.sqrt((nuc_unc_814/nuc_flux_814)**2+(tot_unc_814/tot_flux_814)**2)

# F160W
nuc_flux_160 = flux(data['m160'])
nuc_unc_160 = nuc_flux_160 / 1.086 * data['u160']

tot_flux_160 = flux(dtot['m160'])
tot_unc_160 = tot_flux_160 / 1.086 * np.sqrt((dtot['u160'])**2 + 0.03**2)

frac_160 = nuc_flux_160 / tot_flux_160
func_160 = frac_160 * np.sqrt((nuc_unc_160/nuc_flux_160)**2+(tot_unc_160/tot_flux_160)**2)

for i in range(0,len(data)):
    print(data['Galaxy'][i], frac_475[i], func_475[i])
    print(frac_814[i], func_814[i])
    print(frac_160[i], func_160[i])

name = 'light_mass_frac.pdf'

with PdfPages(name) as pdf:
    fig = plt.figure()

    plt.scatter(frac_475, mass_frac, color='blue')
    #plt.errorbar(frac_475, mass_frac, yerr=[diff_frac_lo, diff_frac_hi], xerr=func_475, color='blue', fmt='o')
    #plt.errorbar(frac_475, mass_frac, xerr=func_814, fmt='o', color='blue')
    plt.scatter(frac_814, mass_frac, color='green')
    #plt.errorbar(frac_814, mass_frac, xerr=func_814, fmt='o', color='green')
    plt.scatter(frac_160, mass_frac, color='red')
    #plt.errorbar(frac_160, mass_frac, xerr=func_160, fmt='o', color='red')
    plt.errorbar(frac_160, mass_frac, yerr=[diff_frac_lo, diff_frac_hi], xerr=func_160, color='red', fmt='o')
    plt.xlabel('Fraction of light associated with compact starburst')
    plt.ylabel('Fraction of stellar mass associated with nuclear starburst')
    plt.xlim(0,1)
    plt.ylim(0,1)
    
    for i in range(0,len(galaxies)):
        plt.text(frac_475[i], mass_frac[i], galaxies[i])
    
    
    pdf.savefig()
    plt.close



# open the pdf file
os.system('open %s &' % name)
