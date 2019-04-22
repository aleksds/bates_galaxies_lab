import numpy as np
import os
from astropy.constants import G
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy import units as u

#constructing cosmological constant
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# total SSP stellar mass
m_ssp_lo = np.array([1.07e10, 6.68e9, 1.06e10, 1.06e10, 7.72e9, 1.07e10, 1.05e10, 1.07e10, 1.06e10, 7.16e9, 1.07e10, 1.45e10]) * const.M_sun #default ssp stellar mass values for 16th percentile
m_ssp_best = np.array([1.24e10, 8.69e9, 1.23e10, 1.24e10, 1.22e10, 1.24e10, 1.23e10, 1.24e10, 1.24e10, 1.18e10, 1.24e10, 1.82e10]) * const.M_sun #50th percentile
m_ssp_up = np.array([1.44e10, 1.15e10, 1.43e10, 1.44e10, 1.43e10, 1.44e10, 1.42e10, 1.43e10, 1.44e10, 1.4e10, 1.43e10, 2.23e10]) * const.M_sun #84th percentile

# total calzetti dust law stellar mass
m_cal_lo = np.array([7.16e9, 8.16e9, 2.13e9, 7.64e9, 9.34e9, 1.07e10, 7.79e9, 2.5e9, 7.95e9, 4.45e9, 2.32e9, 1.27e9]) * const.M_sun #calzetti 16th%
m_cal_best = np.array([1.17e10, 1.45e10, 2.68e9, 1.26e10, 1.47e10, 1.16e10, 1.31e10, 1.22e10, 1.06e10, 3.97e10, 9.29e9, 1.25e10]) * const.M_sun #50th percentile
m_cal_up = np.array([1.37e10, 1.77e10, 1.08e10, 1.66e10, 1.76e10, 1.25e10, 1.67e10, 1.5e10, 1.22e10, 5.45e10, 1.38e10, 1.91e10]) * const.M_sun #calzetti 84th%

# total SFH tau model stellar mass
m_tau_lo = np.array([4.26e9, 4.33e9, 4.28e9, 4.25e9, 4.32e9, 4.34e9, 4.22e9, 4.3e9, 4.79e9, 4.27e9, 4.33e9, 4.33e9]) * const.M_sun #tau 16th%
m_tau_best = np.array([5.56e9, 5.92e9, 5.79e9, 5.53e9, 5.73e9, 6.04e9, 5.38e9, 5.66e9, 6.41e9, 5.6e9, 5.8e9, 5.76e9]) * const.M_sun #50th percentile
m_tau_up = np.array([9.11e9, 9.26e9, 9.26e9, 8.9e9, 9.08e9, 9.36e9, 8.89e9, 9.14e9, 9.51e9, 9.16e9, 9.3e9, 9.23e9]) * const.M_sun #tau 84th%

# total delayed tau model stellar mass
m_dtau_lo = np.array([4.99e9, 5.03e9, 5.02e9, 5.07e9, 5.04e9, 5.08e9, 5.07e9, 5.06e9, 5.49e9, 5.04e9, 5.48e9, 7.58e9]) * const.M_sun #dtau 16th%
m_dtau_best = np.array([6.04e9, 6.04e9, 6.08e9, 6.12e9, 6.07e9, 6.14e9, 6.15e9, 6.07e9, 6.6e9, 6.09e9, 6.6e9, 9.79e9]) * const.M_sun #50th percentile
m_dtau_up = np.array([7.73e8, 7.70e9, 7.73e9, 7.94e9, 7.74e9, 8e9, 8.01e9, 7.74e9, 8.31e9, 7.73e9, 8.26e9, 1.21e10]) * const.M_sun #dtau 84th%

#constructing effective radius limits, redshift, and galaxy notation
values = ascii.read('../autogalfit/ad_mag_size_table.dat')
re_hi_arc = np.array(values['re_large']) * 0.025 * u.arcsec
re_lo_arc = np.array(values['re_small']) * 0.025 * u.arcsec
re_best_arc = np.array(values['re']) * 0.025 * u.arcsec
z = np.array([0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752])
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']

vesc_ssp_best = np.zeros(len(galaxies))
vesc_ssp_lo = np.zeros(len(galaxies))
vesc_ssp_up = np.zeros(len(galaxies))
vesc_cal_up = np.zeros(len(galaxies))
vesc_cal_lo = np.zeros(len(galaxies))
vesc_cal_best = np.zeros(len(galaxies))
vesc_tau_best = np.zeros(len(galaxies))
vesc_tau_up = np.zeros(len(galaxies))
vesc_tau_lo = np.zeros(len(galaxies))
vesc_dtau_up = np.zeros(len(galaxies))
vesc_dtau_lo = np.zeros(len(galaxies))
vesc_dtau_best = np.zeros(len(galaxies))

#finding escape velocity limits for each galaxy
for i in range(0, len(z)):
    re_hi_kpc = re_hi_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])
    re_lo_kpc = re_lo_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])
    re_best_kpc = re_best_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])

    vesc_ssp_lo[i] = np.sqrt(G * m_ssp_lo[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_ssp_best[i] = np.sqrt(G * m_ssp_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_ssp_up[i] = np.sqrt(G * m_ssp_up[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_cal_lo[i] = np.sqrt(G * m_cal_lo[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_cal_best[i] = np.sqrt(G * m_cal_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_cal_up[i] = np.sqrt(G * m_cal_up[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_tau_lo[i] = np.sqrt(G * m_tau_lo[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_tau_best[i] = np.sqrt(G * m_tau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_tau_up[i] = np.sqrt(G * m_tau_up[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_dtau_lo[i] = np.sqrt(G * m_dtau_lo[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_dtau_best[i] = np.sqrt(G * m_dtau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_dtau_up[i] = np.sqrt(G * m_dtau_up[i] / re_lo_kpc).to('km/s') * u.s / u.km
    print(galaxies[i], vesc_ssp_lo,vesc_ssp_best,vesc_ssp_up,vesc_cal_lo,vesc_cal_best,
          vesc_cal_up,vesc_tau_lo,vesc_tau_best,vesc_tau_up,vesc_dtau_lo, vesc_dtau_best,vesc_dtau_up)

#plotting escape velocities

#outflow velocities
vflow = np.array([1228,1206,2470,1778,1828,1830,875,1211,829,2416,1456,606])

#ssp plot
filename = 'kv_vflow_ssp.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_ssp_best, s=30, marker='o', c='blue')
    lower_error_ssp = vesc_ssp_lo
    upper_error_ssp = vesc_ssp_up
    plt.errorbar(vflow, vesc_ssp_best, yerr=[lower_error_ssp, upper_error_ssp], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#cal plot
filename = 'kv_vflow_cal.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_cal_best, s=30, marker='o', c='blue')
    lower_error_cal = vesc_cal_lo
    upper_error_cal = vesc_cal_up
    plt.errorbar(vflow, vesc_cal_best, yerr=[lower_error_cal, upper_error_cal], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#tau plot
filename = 'kv_vflow_tau.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_tau_best, s=30, marker='o', c='blue')
    lower_error_tau = vesc_tau_lo
    upper_error_tau = vesc_tau_up
    plt.errorbar(vflow, vesc_tau_best, yerr=[lower_error_tau, upper_error_tau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#dtau plot
filename = 'kv_vflow_dtau.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_dtau_best, s=30, marker='o', c='blue')
    lower_error_dtau = vesc_dtau_lo
    upper_error_dtau = vesc_dtau_up
    plt.errorbar(vflow, vesc_dtau_best, yerr=[lower_error_dtau, upper_error_dtau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#setting up values for error bars using only stellar mass uncertainties

vesc_ssp_best = np.zeros(len(galaxies))
vesc_ssp_lo = np.zeros(len(galaxies))
vesc_ssp_up = np.zeros(len(galaxies))
vesc_cal_up = np.zeros(len(galaxies))
vesc_cal_lo = np.zeros(len(galaxies))
vesc_cal_best = np.zeros(len(galaxies))
vesc_tau_best = np.zeros(len(galaxies))
vesc_tau_up = np.zeros(len(galaxies))
vesc_tau_lo = np.zeros(len(galaxies))
vesc_dtau_up = np.zeros(len(galaxies))
vesc_dtau_lo = np.zeros(len(galaxies))
vesc_dtau_best = np.zeros(len(galaxies))

for i in range(0, len(z)):
    re_best_kpc = re_best_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])

    vesc_ssp_lo[i] = np.sqrt(G * m_ssp_lo[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_ssp_best[i] = np.sqrt(G * m_ssp_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_ssp_up[i] = np.sqrt(G * m_ssp_up[i] / re_best_kpc).to('km/s') * u.s / u.km

    vesc_cal_lo[i] = np.sqrt(G * m_cal_lo[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_cal_best[i] = np.sqrt(G * m_cal_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_cal_up[i] = np.sqrt(G * m_cal_up[i] / re_best_kpc).to('km/s') * u.s / u.km

    vesc_tau_lo[i] = np.sqrt(G * m_tau_lo[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_tau_best[i] = np.sqrt(G * m_tau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_tau_up[i] = np.sqrt(G * m_tau_up[i] / re_best_kpc).to('km/s') * u.s / u.km

    vesc_dtau_lo[i] = np.sqrt(G * m_dtau_lo[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_dtau_best[i] = np.sqrt(G * m_dtau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_dtau_up[i] = np.sqrt(G * m_dtau_up[i] / re_best_kpc).to('km/s') * u.s / u.km
    print(galaxies[i], vesc_ssp_lo,vesc_ssp_best,vesc_ssp_up,vesc_cal_lo,vesc_cal_best,
          vesc_cal_up,vesc_tau_lo,vesc_tau_best,vesc_tau_up,vesc_dtau_lo, vesc_dtau_best,vesc_dtau_up)

#ssp plot with radius constant
filename = 'kv_vflow_ssp_mass_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_ssp_best, s=30, marker='o', c='blue')
    lower_error_ssp = vesc_ssp_lo
    upper_error_ssp = vesc_ssp_up
    plt.errorbar(vflow, vesc_ssp_best, yerr=[lower_error_ssp, upper_error_ssp], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#cal plot with radius constant
filename = 'kv_vflow_cal_mass_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_cal_best, s=30, marker='o', c='blue')
    lower_error_cal = vesc_cal_lo
    upper_error_cal = vesc_cal_up
    plt.errorbar(vflow, vesc_cal_best, yerr=[lower_error_cal, upper_error_cal], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#tau plot with radius constant
filename = 'kv_vflow_tau_mass_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_tau_best, s=30, marker='o', c='blue')
    lower_error_tau = vesc_tau_lo
    upper_error_tau = vesc_tau_up
    plt.errorbar(vflow, vesc_tau_best, yerr=[lower_error_tau, upper_error_tau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#dtau plot with radius constant
filename = 'kv_vflow_dtau_mass_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_dtau_best, s=30, marker='o', c='blue')
    lower_error_dtau = vesc_dtau_lo
    upper_error_dtau = vesc_dtau_up
    plt.errorbar(vflow, vesc_dtau_best, yerr=[lower_error_dtau, upper_error_dtau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#setting up values for error bars using only effective radius uncertainties

vesc_ssp_best = np.zeros(len(galaxies))
vesc_ssp_lo = np.zeros(len(galaxies))
vesc_ssp_up = np.zeros(len(galaxies))
vesc_cal_up = np.zeros(len(galaxies))
vesc_cal_lo = np.zeros(len(galaxies))
vesc_cal_best = np.zeros(len(galaxies))
vesc_tau_best = np.zeros(len(galaxies))
vesc_tau_up = np.zeros(len(galaxies))
vesc_tau_lo = np.zeros(len(galaxies))
vesc_dtau_up = np.zeros(len(galaxies))
vesc_dtau_lo = np.zeros(len(galaxies))
vesc_dtau_best = np.zeros(len(galaxies))

for i in range(0, len(z)):
    re_hi_kpc = re_hi_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])
    re_lo_kpc = re_lo_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])
    re_best_kpc = re_best_arc[i] / cosmo.arcsec_per_kpc_proper(z[i])

    vesc_ssp_lo[i] = np.sqrt(G * m_ssp_best[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_ssp_best[i] = np.sqrt(G * m_ssp_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_ssp_up[i] = np.sqrt(G * m_ssp_best[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_cal_lo[i] = np.sqrt(G * m_cal_best[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_cal_best[i] = np.sqrt(G * m_cal_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_cal_up[i] = np.sqrt(G * m_cal_best[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_tau_lo[i] = np.sqrt(G * m_tau_best[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_tau_best[i] = np.sqrt(G * m_tau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_tau_up[i] = np.sqrt(G * m_tau_best[i] / re_lo_kpc).to('km/s') * u.s / u.km

    vesc_dtau_lo[i] = np.sqrt(G * m_dtau_best[i] / re_hi_kpc).to('km/s') * u.s / u.km
    vesc_dtau_best[i] = np.sqrt(G * m_dtau_best[i] / re_best_kpc).to('km/s') * u.s / u.km
    vesc_dtau_up[i] = np.sqrt(G * m_dtau_best[i] / re_lo_kpc).to('km/s') * u.s / u.km
    print(galaxies[i], vesc_ssp_lo,vesc_ssp_best,vesc_ssp_up,vesc_cal_lo,vesc_cal_best,
          vesc_cal_up,vesc_tau_lo,vesc_tau_best,vesc_tau_up,vesc_dtau_lo, vesc_dtau_best,vesc_dtau_up)

#ssp plot with stellar mass constant
filename = 'kv_vflow_ssp_re_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_ssp_best, s=30, marker='o', c='blue')
    lower_error_ssp = vesc_ssp_lo
    upper_error_ssp = vesc_ssp_up
    plt.errorbar(vflow, vesc_ssp_best, yerr=[lower_error_ssp, upper_error_ssp], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#cal plot with stellar mass constant
filename = 'kv_vflow_cal_re_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_cal_best, s=30, marker='o', c='blue')
    lower_error_cal = vesc_cal_lo
    upper_error_cal = vesc_cal_up
    plt.errorbar(vflow, vesc_cal_best, yerr=[lower_error_cal, upper_error_cal], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#tau plot with stellar mass constant
filename = 'kv_vflow_tau_re_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_tau_best, s=30, marker='o', c='blue')
    lower_error_tau = vesc_tau_lo
    upper_error_tau = vesc_tau_up
    plt.errorbar(vflow, vesc_tau_best, yerr=[lower_error_tau, upper_error_tau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)

#dtau plot with stellar mass constant
filename = 'kv_vflow_dtau_re_unc.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(vflow, vesc_dtau_best, s=30, marker='o', c='blue')
    lower_error_dtau = vesc_dtau_lo
    upper_error_dtau = vesc_dtau_up
    plt.errorbar(vflow, vesc_dtau_best, yerr=[lower_error_dtau, upper_error_dtau], fmt='o')
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Estimated Escape Velocity')
    plt.xlim(0,3000)
    plt.ylim(0,3000)
    x = np.arange(0,3001)
    print (x)
    y = x
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)