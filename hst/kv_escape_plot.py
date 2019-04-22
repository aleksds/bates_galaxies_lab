#calculates escape velocity from stellar masses obtained from prospector and radii obtained from galfit. does so for 475 and 814 using the best value prosepector gave us for stellar mass and the two percentile values
import numpy as np
import os
from astropy.constants import G
from astropy import constants as const
from astropy import units as u
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# stellar mass
masslogsol = np.array([10.48,10.58,10.12,10.61,10.54,11.09,10.29,10.23,10.98,11.35,10.67,10.11]) #in standard galactic order
mass = 10**masslogsol * const.M_sun

#getting radius in meters from pc
rem = np.zeros([12,2])
repc = np.array([[84,134],[175,149],[86,83],[158,74],[62,179],[164,239],[170,81],[77,108],[215,1032],[560,886],[98,178],[91,189]])
re = repc * const.pc

#now plugging into escape velocity formula to get v_esc
vesc_475 = np.sqrt(2*G*mass/re[:,0]).to('km/s')
vesc_814 = np.sqrt(2*G*mass/re[:,1]).to('km/s')

#now getting escape velocity from the 16th and 18th percentile values for the stellar mass from prospector

#need to get the masses into kgs
#first for the 16th percentile
masslogsolsix = np.array([10.29,10.40,9.88,10.41,10.36,10.92,10.10,10.00,10.79,11.18,10.48,9.89])
mass_16 = 10**masslogsolsix * const.M_sun
vesc_16_475 = np.sqrt(2*G*mass_16/re[:,0]).to('km/s')
vesc_16_814 = np.sqrt(2*G*mass_16/re[:,1]).to('km/s')


#now for the 84th percentile
masslogsoleight = np.array([10.70,10.75,10.32,10.79,10.73,11.27,10.48,10.42,11.18,11.51,10.91,10.38])
mass_84 = 10**masslogsoleight * const.M_sun
vesc_84_475 = np.sqrt(2*G*mass_84/re[:,0]).to('km/s')
vesc_84_814 = np.sqrt(2*G*mass_84/re[:,1]).to('km/s')

vflow = np.array([1228,1206,2470,1778,1828,1830,875,1211,829,2416,1456,606])

filename = 'kv_vflowvsescvel.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_475, s=30, marker='v', c='blue', label='F475W')
    lower_error_475 = (vesc_475 - vesc_16_475) * u.s / u.km
    upper_error_475 = (vesc_84_475 - vesc_475) * u.s / u.km
    plt.errorbar(vflow, vesc_475 * u.s / u.km, yerr=[lower_error_475, upper_error_475], fmt='v')
    plt.scatter(vflow, vesc_814, s=30, marker='o', c='green', label = 'F814W')
    lower_error_814 = (vesc_814 - vesc_16_814) * u.s / u.km
    upper_error_814 = (vesc_84_814 - vesc_814) * u.s / u.km
    plt.errorbar(vflow, vesc_814 * u.s / u.km, yerr=[lower_error_814, upper_error_814], fmt='o', c='green')
    plt.legend(loc='upper left', prop={'size': 12})
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)


jc_values = ascii.read('autogalfit/jc_mag_size_table.dat')
re_pix = np.array(jc_values['re'])
z = np.array([0.603,0.459,0.712,0.514,0.467,0.451,0.658,0.608,0.402,0.449,0.728,0.752])
re_kpc = np.zeros(len(z))

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
for i in range(0,len(z)):
    re_kpc[i] = re_pix[i] / cosmo.arcsec_per_kpc_proper(z[i]) * 0.025 * u.arcsec / u.kpc
    print(re_pix[i], re_kpc[i])

re_kpc = re_kpc * u.kpc
mass_test = np.zeros(12) + 1e10*const.M_sun
vesc_test = np.sqrt(2*G*mass_test/re_kpc).to('km/s')

filename = 'ad_vflow_test.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_test, s=30, marker='o', c='blue')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)


mass_50_new = np.array([10.23,10.39,10.18,10.34,10.05,10.61,10.08,10.01,10.34,10.91,10.21,10.37])
mass_84_diff = np.array([0.12, 0.19, 0.17, 0.20, 0.19, 0.15, 0.16, 0.18, 0.19, 0.15, 0.18, 0.16])
mass_16_diff = np.array([0.15, 0.18, 0.14, 0.14, 0.15, 0.17, 0.14, 0.12, 0.21, 0.17, 0.12, 0.20])

mass_new = 10.**(mass_50_new)*const.M_sun
mass_84_new = 10.**(mass_50_new + mass_84_diff)*const.M_sun
mass_16_new = 10.**(mass_50_new - mass_16_diff)*const.M_sun

vesc_new = np.sqrt(2*G*mass_new/re_kpc).to('km/s')
vesc_84_new = np.sqrt(2*G*mass_84_new/re_kpc).to('km/s')
vesc_16_new = np.sqrt(2*G*mass_16_new/re_kpc).to('km/s')

filename = 'ad_vflow_new.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_new, s=30, marker='o', c='blue')
    lower_error_new = (vesc_new - vesc_16_new) * u.s / u.km
    upper_error_new = (vesc_84_new - vesc_new) * u.s / u.km
    plt.errorbar(vflow, vesc_new * u.s / u.km, yerr=[lower_error_new, upper_error_new], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

mass_50_ssp =  np.array([9.38,9.26,9.73,9.36,9.30,9.41,9.34,9.54,9.36,9.52,9.30,9.71])
mass_84_diff = np.array([0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20])
mass_16_diff = np.array([0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20,0.20])

mass_ssp = 10.**(mass_50_ssp)*const.M_sun
mass_84_ssp = 10.**(mass_50_ssp + mass_84_diff)*const.M_sun
mass_16_ssp = 10.**(mass_50_ssp - mass_16_diff)*const.M_sun

vesc_ssp = np.sqrt(2*G*mass_ssp/re_kpc).to('km/s')
vesc_84_ssp = np.sqrt(2*G*mass_84_ssp/re_kpc).to('km/s')
vesc_16_ssp = np.sqrt(2*G*mass_16_ssp/re_kpc).to('km/s')

filename = 'ad_vflow_ssp.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_ssp, s=30, marker='o', c='blue')
    lower_error_ssp = (vesc_ssp - vesc_16_ssp) * u.s / u.km
    upper_error_ssp = (vesc_84_ssp - vesc_ssp) * u.s / u.km
    plt.errorbar(vflow, vesc_ssp * u.s / u.km, yerr=[lower_error_ssp, upper_error_ssp], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

r = 100 * const.pc
m = 2e9 * const.M_sun
print(np.sqrt(2.*G*m/r).to('km/s'))
