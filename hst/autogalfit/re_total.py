from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u


# function to return a flux in maggies given an AB magnitude
def flux(mag):
    flux = 10. ** (mag / (-2.5))
    return flux

data = ascii.read('bgl_phot.dat')
dtot = ascii.read('bgl_tot.dat')

galaxies =        ['J0826','J0901','J0905','J0944','J1107','J1219','J1341','J1506','J1558','J1613','J2116','J2140']

re_best = np.array([0.0151, 0.0149, 0.0105, 0.0099, 0.0156, 0.0257, 0.0127, 0.0118, 0.0387, 0.1289, 0.0216, 0.0145]) 
re_unc = np.array([0.0031, 0.0033, 0.0027, 0.0030, 0.0041, 0.0038, 0.0023, 0.0025, 0.0064, 0.0157, 0.0046, 0.0045])
re_best_arc = re_best * u.arcsec

z = np.array([0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752])

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

re_best_kpc = re_best_arc / cosmo.arcsec_per_kpc_proper(z)
re_best_kpc_hi = (re_best_arc + re_unc * u.arcsec) / cosmo.arcsec_per_kpc_proper(z)
re_best_kpc_lo = (re_best_arc - re_unc * u.arcsec) / cosmo.arcsec_per_kpc_proper(z)



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

# hack for J1558
frac_814[8] = 0.6

# F160W
nuc_flux_160 = flux(data['m160'])
nuc_unc_160 = nuc_flux_160 / 1.086 * data['u160']

tot_flux_160 = flux(dtot['m160'])
tot_unc_160 = tot_flux_160 / 1.086 * np.sqrt((dtot['u160'])**2 + 0.03**2)

frac_160 = nuc_flux_160 / tot_flux_160
func_160 = frac_160 * np.sqrt((nuc_unc_160/nuc_flux_160)**2+(tot_unc_160/tot_flux_160)**2)

# define de Vaucouleurs profile
r = (np.arange(10000)+1) / 100.
re = 1.
surf = np.exp(-7.669*((r/re)**0.25-1))
flux = np.zeros(len(r))
flux[0] = surf[0] * 2 * np.pi * r[0] * r[0]

for i in range(1,len(r)):
    flux[i] = flux[i-1] + surf[i] * 2. * np.pi * r[i] * (r[i] - r[i-1])

flux = flux / np.max(flux)

extrap = np.zeros(len(data))
extrap_min = np.zeros(len(data))
extrap_max = np.zeros(len(data))

fac_ext = np.zeros(len(data))
fac_ext_min = np.zeros(len(data))
fac_ext_max = np.zeros(len(data))


for i in range(0,len(data)):
    print(data['Galaxy'][i])
    print('F475W:', frac_475[i], func_475[i])
    print('F814W:', frac_814[i], func_814[i])
    print('F160W:', frac_160[i], func_160[i])

    extrap[i] = 0.5/frac_814[i]
    dif = np.min(np.abs(extrap[i]-flux))
    good = np.where(np.abs(extrap[i]-flux) == dif)
    fac_ext[i] = r[good]

    print('need to extraopolate from 0.5 to:', extrap[i])
    print(flux[good], fac_ext[i])

    extrap_min[i] = 0.5/(frac_814[i]+func_814[i])
    dif = np.min(np.abs(extrap_min[i]-flux))
    good = np.where(np.abs(extrap_min[i]-flux) == dif)
    fac_ext_min[i] = r[good]
    
    print('at minimum, need to extrapolate to:', extrap_min[i])
    print(flux[good], fac_ext_min[i])
    
    extrap_max[i] = 0.5/(frac_814[i]-func_814[i])
    dif = np.min(np.abs(extrap_max[i]-flux))
    good = np.where(np.abs(extrap_max[i]-flux) == dif)
    fac_ext_max[i] = r[good]

    print('at maximum, need to extrapolate to:', extrap_max[i])
    print(flux[good], fac_ext_max[i])

 
re_best_kpc_cor = re_best_kpc * fac_ext
re_best_kpc_cor_hi = re_best_kpc_hi * fac_ext_max
re_best_kpc_cor_lo = re_best_kpc_lo * fac_ext_min




