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

#ext_best_mass =   np.array([10.66, 10.55,  10.12,  10.54,  10.67,  10.99,  10.43,  11.03,  10.63,  10.65,  10.98,  10.42])
#ext_lo_mass =     np.array([0.32,  0.26,   0.18,   0.14,   0.17,   0.16,   0.44,   0.16,   0.15,   0.47,   0.27,   0.25])
#ext_up_mass =     np.array([0.21,  0.19,   0.17,   0.13,   0.13,   0.10,   0.43,   0.15,   0.13,   0.47,   0.19,   0.20])

#                        10.8>11.0        >10.95 10.5>10.7                     10.7>10.9 10.6>10.8 11.2>11.4      10.9>11.1
#                           JM      JM      JM      JM      JM     JM        JM     CT      JM      JM     JM      JM
tremonti_mass =  np.array([11.03, 10.98,  10.74,  10.69,  10.73,  10.71,  10.98,  10.96,  10.91,  10.98,  10.92,  11.18])
tremonti_vavg = np.array([1131.2,1210.1, 2395.3, 1172.4, 1359.8, 1543.7,  766.9, 1378.9,  807.1, 1437.5, 1056.2,  343.8])   
tremonti_vmax = np.array([1456.2,1575.4, 2884.8, 1878.3, 2015.6, 1962.6, 1995.2, 1768.8, 1151.2, 2374.1, 1915.1,  950.4])   

moustakas_mass = np.array([11.31, 11.60,  11.27,  11.15,  11.01,  11.30,  10.90,  11.62,  11.17,  11.44,  11.10,  11.11])

tot_best_mass =  np.array([10.71, 10.59,  10.75,  10.74,  10.68,  11.34,  10.61,  10.65,  10.97,  11.34,  11.16,  11.12])
tot_lo_mass =     np.array([0.01,  0.01,   0.02,   0.03,   0.02,   0.02,   0.02,   0.03,   0.02,   0.01,   0.05,   0.05])
tot_up_mass =     np.array([0.02,  0.02,   0.02,   0.04,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.04,   0.05])

#                           0.92   0.95    0.10    0.83    0.94    0.65    0.13    0.73    0.20    0.17    0.76    0.86
new_best_mass =  np.array([10.60, 10.50,  10.72,  10.58,  10.60,  10.94,  10.76,  10.65,  10.65,  11.04,  10.92,  10.95])
new_lo_mass =     np.array([0.01,  0.01,   0.01,   0.02,   0.01,   0.03,   0.01,   0.03,   0.02,   0.02,   0.04,   0.02])
new_up_mass =     np.array([0.02,  0.01,   0.01,   0.03,   0.02,   0.05,   0.01,   0.03,   0.03,   0.02,   0.06,   0.02])

#                           0.91,  0.93,   0.93    0.85    0.93    0.75    0.96    0.78    0.40,   0.18,   0.77,   0.88
you_best_mass =  np.array([10.62, 10.51,  10.68,  10.60,  10.62,  10.97,  10.57,  10.67,  10.72,  11.14,  10.99,  10.96])
you_lo_mass   =  np.array([ 0.03,  0.02,   0.02,   0.05,   0.03,   0.06,   0.01,   0.06,   0.09,   0.05,   0.09,   0.03])
you_up_mass   =  np.array([ 0.04,  0.03,   0.04,   0.06,   0.03,   0.07,   0.02,   0.08,   0.05,   0.06,   0.09,   0.06])

#                           0.34   0.44    0.10    0.40    0.39    0.30    0.12    0.39    0.16    0.24,   0.33    0.12
#                           0.39   0.35    0.49    0.34    0.32    0.28    0.58    0.30    0.15    0.37    0.51    0.39   
lew_best_mass =  np.array([10.68, 10.56,  10.72,  10.67,  10.60,  11.00,  10.75,  10.71,  10.66,  11.08,  10.97,  10.97]) 
lew_lo_mass =    np.array([ 0.03,  0.01,   0.02,   0.03,   0.03,   0.03,   0.01,   0.03,   0.03,   0.05,   0.05,   0.03])
lew_up_mass =    np.array([ 0.02,  0.02,   0.02,   0.04,   0.06,   0.05,   0.01,   0.03,   0.03,   0.04,   0.07,   0.03])

#                           0.45   0.45    0.45    0.44    0.47     0.40   0.47    0.43,   0.31,   0.10,   0.41,   0.45 
#                           0.23   0.19    0.27    0.25    0.17     0.22   0.44    0.22    0.12    0.30    0.39    0.24
dew_best_mass =  np.array([10.90, 10.81,  10.98,  10.80,  10.89,   11.11, 10.86,  10.84,  10.77,  11.13,  11.11,  11.16])
dew_lo_mass =    np.array([ 0.03,  0.03,   0.03,   0.05,   0.04,    0.05,  0.02,   0.04,   0.05,   0.04,   0.07,   0.06])
dew_up_mass =    np.array([ 0.06,  0.05,   0.05,   0.06,   0.04,    0.06,  0.04,   0.05,   0.06,   0.05,   0.09,   0.05])

#                           0.46   0.45    0.09     0.44   0.10
#                           0.35   0.19    0.37     0.31   0.06
#cew_best_mass =  np.array([10.73, 10.82,  10.84,  10.71  11.32
#cew_lo_mass =    np.array([ 0.03,  0.04,   0.02,   0.03   0.12
#cew_hi_mass =    np.array([ 0.03,  0.05,   0.02,   0.03   0.10

#                           0.36   0.36    0.37    0.36   0.12     0.34    0.38    0.35   0.28    0.17    0.33    0.37
#                           0.19   0.16    0.22    0.21   0.08     0.20    0.36    0.20   0.11    0.32    0.35    0.19
# 11.22 vs 10.87 for J1107
few_best_mass =  np.array([10.99, 10.90,  11.06,  10.88, 11.22,   11.15,  10.95,  10.88, 10.79,  11.14,  11.14,  11.28])
few_lo_mass =    np.array([ 0.03,  0.03,   0.03,   0.08,  0.07,    0.06,   0.02,   0.04,  0.04,   0.04,   0.06,   0.08])
few_up_mass =    np.array([ 0.05,  0.06,   0.05,   0.08,  0.07,    0.07,   0.03,   0.06,  0.08,   0.06,   0.08,   0.06])


#tot_lo_mass = np.sqrt(tot_lo_mass**2 + 0.1**2)
#tot_up_mass = np.sqrt(tot_up_mass**2 + 0.1**2)

mass_frac = 10**(nuc_best_mass) / (10**(tot_best_mass))
mass_frac_hi = 10**(nuc_best_mass+nuc_up_mass) / (10**(tot_best_mass-tot_lo_mass))
mass_frac_lo = 10**(nuc_best_mass-nuc_lo_mass) / (10**(tot_best_mass+tot_up_mass))

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

name = 'light_tot_frac.pdf'


numbers = np.arange(110)/100.

with PdfPages(name) as pdf:
    #fig = plt.figure()
    #
    #plt.scatter(frac_475, mass_frac, color='blue', marker='*', label='F475W')
    ##plt.errorbar(frac_475, mass_frac, yerr=[diff_frac_lo, diff_frac_hi], xerr=func_475, color='blue', fmt='*', elinewidth=1)
    #plt.errorbar(frac_475, mass_frac, xerr=func_814, fmt='*', color='blue',elinewidth=1)
    #plt.scatter(frac_814, mass_frac, color='green', marker='^', label='F814W')
    ##plt.errorbar(frac_814, mass_frac, xerr=func_814, fmt='^', color='green',elinewidth=1)
    #plt.errorbar(frac_814, mass_frac, yerr=[diff_frac_lo, diff_frac_hi], xerr=func_814, color='green', fmt='^', elinewidth=1)
    #plt.scatter(frac_160, mass_frac, color='red', marker='o', label='F160W')
    #plt.errorbar(frac_160, mass_frac, xerr=func_160, fmt='o', color='red', elinewidth=1)
    ##plt.errorbar(frac_160, mass_frac, yerr=[diff_frac_lo, diff_frac_hi], xerr=func_160, color='red', fmt='o', elinewidth=1)
    #plt.plot(numbers, numbers, color='black', label='mass fraction = light fraction')
    #plt.xlabel('Fraction of light associated with compact starburst')
    #plt.ylabel('Fraction of stellar mass associated with compact starburst')
    #plt.xlim(0.35,1)
    #plt.ylim(0,1)
    #
    #
    #for i in range(0,len(galaxies)):
    #    plt.text(frac_475[i], mass_frac[i], galaxies[i])
    #    plt.plot(np.array([frac_475[i], frac_814[i], frac_160[i]]), np.zeros(3)+mass_frac[i], color='gray', linestyle=':')
    #
    #
    #plt.plot(np.array([frac_475[0], frac_814[0], frac_160[0]]), np.zeros(3)+mass_frac[0], color='gray', linestyle=':', label='range of light fractions')
    #plt.legend()
    #    
    #pdf.savefig()
    #plt.close()
    #
    #fig = plt.figure()
    #
    #
    #plt.scatter(tot_best_mass, mass_frac)
    #
    #plt.errorbar(tot_best_mass, mass_frac, xerr=np.zeros(12)+0.1, yerr=[diff_frac_lo, diff_frac_hi], elinewidth=1, fmt='o')
    #
    #plt.xlim(10.5,11.5)
    #plt.ylim(0,1)
    #
    #pdf.savefig()
    #plt.close()


    fig = plt.figure()

    yo = np.arange(101)/50+9.5
    
    #plt.scatter(tot_best_mass, nuc_best_mass)
    #plt.scatter(tremonti_mass, nuc_best_mass)
    plt.scatter(dew_best_mass, nuc_best_mass)
    #plt.scatter(few_best_mass, nuc_best_mass)
    #plt.errorbar(tot_best_mass, nuc_best_mass, yerr=[np.sqrt(nuc_lo_mass**2+0.1**2), np.sqrt(nuc_up_mass**2+0.1**2)], xerr=[np.sqrt(tot_lo_mass**2+0.1**2), np.sqrt(tot_up_mass**2+0.1**2)], fmt='o', elinewidth=1)
    #plt.errorbar(tremonti_mass, nuc_best_mass, yerr=[np.sqrt(nuc_lo_mass**2+0.1**2), np.sqrt(nuc_up_mass**2+0.1**2)], xerr=[np.sqrt(tot_lo_mass**2+0.1**2), np.sqrt(tot_up_mass**2+0.1**2)], fmt='o', elinewidth=1)
    plt.errorbar(dew_best_mass, nuc_best_mass, yerr=[np.sqrt(nuc_lo_mass**2+0.1**2), np.sqrt(nuc_up_mass**2+0.1**2)], xerr=[np.sqrt(dew_lo_mass**2+0.1**2), np.sqrt(dew_up_mass**2+0.1**2)], fmt='o', elinewidth=1, color='#ff7f0e') # 20200518
    #plt.errorbar(few_best_mass, nuc_best_mass, yerr=[np.sqrt(nuc_lo_mass**2+0.1**2), np.sqrt(nuc_up_mass**2+0.1**2)], xerr=[np.sqrt(few_lo_mass**2+0.1**2), np.sqrt(few_up_mass**2+0.1**2)], fmt='o', elinewidth=1)

    #plt.plot(yo,yo) # 20200518
    #plt.plot(yo,yo, label=r'${M}_{*,central}=M_{*,total}$' )
    plt.plot(yo,yo, label=r'$\mathcal{M}_{*,central}=\mathcal{M}_{*,total}$', linestyle='dashed', color='#2ca02c' ) # 20200518
    #plt.plot(yo,yo+np.log10(1/3), linestyle=':', label=r'$\mathcal{M}_{*,central}=\frac{1}{3}\mathcal{M}_{*,total}$')
    #plt.plot(yo,yo+np.log10(0.20), linestyle=':', label=r'$\mathcal{M}_{*,central}=0.2\mathcal{M}_{*,total}$')
    plt.plot(yo,yo+np.log10(0.25), linestyle=':', label=r'$\mathcal{M}_{*,central}=0.25\mathcal{M}_{*,total}$', color='#1f77b4') # 20200518
    #plt.plot(yo,yo+np.log10(0.35), linestyle='--', label=r'${M}_{*,central}=10\%\ M_{*,total}$')
    #plt.plot(yo,yo+np.log10(0.25), linestyle='-.', label=r'${M}_{*,central}=50\%\ M_{*,total}$')

    
    plt.xlabel(r'$\log_{10}(\mathcal{M}_{*, total\ }/\mathcal{M}_{\odot})$', fontsize=14)
    #ax1.set_xlabel(r'$\log_{10}(\mathcal{M}/\mathcal{M}_{\odot})$')
    #plt.ylabel(r'$\log{(M_{*, central})}$', fontsize=14)
    plt.ylabel(r'$\log_{10}(\mathcal{M}_{*, central\ }/\mathcal{M}_{\odot})$', fontsize=14)

    plt.legend(loc='lower right', fontsize=14)

    #for i in range(0,len(galaxies)):
        #plt.text(tot_best_mass[i], nuc_best_mass[i], galaxies[i])
        #plt.text(tremonti_mass[i], nuc_best_mass[i], galaxies[i])
        #plt.text(dew_best_mass[i], nuc_best_mass[i], galaxies[i])
        #plt.text(few_best_mass[i], nuc_best_mass[i], galaxies[i])

    

    #plt.xlim(10.4,11.3)
    plt.xlim(10.6,11.4)
    plt.ylim(9.6,10.9)
    
    pdf.savefig()
    plt.close()

# open the pdf file
os.system('open %s &' % name)
