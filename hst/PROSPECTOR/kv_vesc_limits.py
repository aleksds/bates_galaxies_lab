import numpy as np
import os
from astropy.constants import G
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy import units as u
import xlrd

#constructing cosmological constant
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

#reading into excel file for stellar mass values
loc = ('../PROSPECTOR/Results worksheet for stellar mass, age, and dust.xlsx')
wb = xlrd.open_workbook(loc)
mass_sheet = wb.sheet_by_index(0)

#setting up stellar mass array
m_ssp_lo = np.zeros([])
m_ssp_best = np.zeros([])
m_ssp_up = np.zeros([])

m_cal_lo = np.zeros([])
m_cal_best = np.zeros([])
m_cal_up = np.zeros([])

m_tau_lo = np.zeros([])
m_tau_best = np.zeros([])
m_tau_up = np.zeros([])

m_dtau_lo = np.zeros([])
m_dtau_best = np.zeros([])
m_dtau_up = np.zeros([])

for i in range(0, mass_sheet.nrows-1):
    # total SSP stellar mass
    m_ssp_lo = np.array([mass_sheet.cell_value(i+1,9)]) * const.M_sun #default ssp stellar mass values for 16th percentile
    m_ssp_best = np.array([mass_sheet.cell_value(i+1,1)]) * const.M_sun #50th percentile
    m_ssp_up = np.array([mass_sheet.cell_value(i+1,5)]) * const.M_sun #84th percentile

    #calzetti dust law
    m_cal_lo = np.array([mass_sheet.cell_value(i+1,10)]) * const.M_sun
    m_cal_best = np.array([mass_sheet.cell_value(i+1,2)]) * const.M_sun
    m_cal_up = np.array([mass_sheet.cell_value(i+1,6)]) * const.M_sun

    #SFH tau model
    m_tau_lo = np.array([mass_sheet.cell_value(i+1,11)]) * const.M_sun
    m_tau_best = np.array([mass_sheet.cell_value(i+1,3)]) * const.M_sun
    m_tau_up = np.array([mass_sheet.cell_value(i+1,7)]) * const.M_sun

    #SFH Delayed tau model
    m_dtau_lo = np.array([mass_sheet.cell_value(i+1,12)]) * const.M_sun
    m_dtau_best = np.array([mass_sheet.cell_value(i+1,4)]) * const.M_sun
    m_dtau_up = np.array([mass_sheet.cell_value(i+1,8)]) * const.M_sun


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
"""
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

os.system('open %s &' % filename)"""