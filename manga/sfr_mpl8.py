# Rebecca Minsley
# April 29th 2019

import os
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy import units as u
from os.path import isdir, join
import csv
import pandas as pd
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)      # define a standard cosmology
mpl8_dir = os.environ['MANGADIR_MPL8']  # directory should be saved in .bash_profile

# data release pipline
drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')            # read in information from drpall file
drpdata = drp[1].data                                       # drpall
manga_ra = drpdata.field('objra')                           # declination
manga_dec = drpdata.field('objdec')
m_survey = drpdata.field('srvymode') == 'MaNGA dither'      # boolean -> True at MaNGA Data
ba = drpdata.field('nsa_sersic_ba')[m_survey]               # sersic profiles of only MANGA galaxies
# nsa
nsa = fits.open(mpl8_dir + '1-nsa_v1_0_1.fits')             # read in information from NSA catalog
nsa_data = nsa[1].data                                      # nsa catalog
nsa_ra = nsa_data.field('RA')
nsa_dec = nsa_data.field('DEC')
th50 = nsa_data.field('ELPETRO_TH50_R')                     # th50 concentration array
th90 = nsa_data.field('ELPETRO_TH90_R')                     # th90 concentration array


# match the two catalogs on the basis of right ascension and declination
c = SkyCoord(ra=manga_ra*u.degree, dec=manga_dec*u.degree)
catalog = SkyCoord(ra=nsa_ra*u.degree, dec=nsa_dec*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)
th50_manga = th50[idx]
th90_manga = th90[idx]
c_manga = (th90_manga/th50_manga)[m_survey]                 # MaNGA galaxy concentration used to define late-type


# late-type edge on
late = c_manga < 2.6                                        # late-type concentrations (<2.6)
edge = ba[late] < 0.3                                       # Defining edge-on Galaxies
stellar_mass = drpdata.field('nsa_sersic_mass')[m_survey][late][edge]             # MaNGA gal stellar mass
log_stellar_mass = np.log10(stellar_mass)
plateifu = drpdata.field('plateifu')[m_survey][late][edge]                        # MaNGA plateifu in h-2 solar masses
nsa_z = drpdata.field('nsa_zdist')[m_survey][late][edge]                          # multiply by c/H0 for Mpc.



# MPA-JHU data
mpa_gal = fits.open(mpl8_dir + 'gal_info_dr7_v5_2.fit')     # MPA-JHU catalog galaxy data
mpa_galdata = mpa_gal[1].data
mpa_ra = mpa_galdata.field('RA')
mpa_dec = mpa_galdata.field('DEC')

mpa_sfr = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')  # MPA-JHU catalog star-formation
mpa_sfrdata = mpa_sfr[1].data


# match the MPA-JHU and MaNGA catalogs on the basis of RA and DEC
bad = mpa_dec < -90
mpa_dec[bad] = -90
mpa_cat = SkyCoord(ra=mpa_ra*u.degree, dec=mpa_dec*u.degree)
idx_mpa, d2d_mpa, d3d_mpa = c.match_to_catalog_sky(mpa_cat)
mpa_sfr = mpa_sfrdata.field('MEDIAN')[idx_mpa]      # MPA-JHU catalog matched data
mpa_sfr_manga = mpa_sfr[m_survey][late][edge]       # MPA-JHU late-type edge-on MaNGA galaxy sfr



 # Data Access Pipeline
dap = fits.open(mpl8_dir + 'dapall-v2_5_3-2.3.0.fits')      # read in information from dapall file
dapdata = dap[1].data                                       # dapall
dap_sfr = dapdata.field('SFR_TOT')                              # sfr rate h_flux lumin from DAP file


good_plates = plateifu
manga_star = np.zeros(len(good_plates))
manga_starb = np.zeros(len(good_plates))
total = len(good_plates)
print('The total number of late-type, edge-on MaNGA galaxies is ' + str(total))


for i in range(5):
    plate, ifu = good_plates[i].split('-')
    name = mpl8_dir+'HYB10-MILESHC-MILESHC/'+str(plate)+'/'+ str(ifu)+'/manga-'+str(plate)+'-'+str(ifu)+'-MAPS-HYB10-MILESHC-MILESHC.fits.gz'
    dap_name = mpl8_dir+ 'manga-' + str(good_plates[i]) + '-MAPS-HYB10-MILESHC-MILESHC.fits.gz'

    if os.path.isfile(name):
        match = np.where((drpdata.field('plate') == int(plate)) & (drpdata.field('ifudsgn') == str(ifu)))
        hdu_drp = fits.open(name)
        print(str(plate) + '-' + str(ifu))

        emlc = channel_dictionary(hdu_drp, 'EMLINE_SFLUX')                      # dictionary w/ emission-line & spectral-index channel names
        dap_ha_sflux = hdu_drp['EMLINE_SFLUX'].data[emlc['Ha-6564'], :, :]      # h-alpha summed emission flux (10^(-17) erg/s/cm2/spaxel)
        dap_hb_sflux = hdu_drp['EMLINE_SFLUX'].data[emlc['Hb-4862'], :, :]

        ha_flux = np.sum(dap_ha_sflux)
        hb_flux = np.sum(dap_hb_sflux)

        # intrinsic h_alpha/h_beta ratio for Case B recombination is 2.86
        color_excess = 0.934 * math.log((ha_flux/hb_flux)/2.86)
        k_halpha = 2.468

        z = cosmo.luminosity_distance(nsa_z[i]).to(u.cm)     #(nsa_z[i](const.c / cosmo.H0.to(u.m/(u.m * u.s)))).to(u.cm)
        obs_ha_flux = (ha_flux * (10 ** (-0.4 * k_halpha * color_excess)))    # in 10^-17 erg/s/cm2/spaxel
        obs_ha_flux = obs_ha_flux * 10**-17 * u.erg/u.s/(u.cm ** 2)

        lum = 4 * math.pi * z**2 * obs_ha_flux
        log_sfr = math.log(lum.value, 10) - 41.27
        manga_star[i] = log_sfr



        lum_b = (4 * math.pi * z**2 * ha_flux) * (((ha_flux/hb_flux)/2.8)**2.36)
        sfrb = lum_b/(10**41.1)
        manga_starb[i] = math.log(sfrb.value, 10)
        # mpl8_sfr = dap_data['SFR_TOT']

with PdfPages(filename) as pdf:

    # stellar mass vs star-formation
    fig = plt.figure()
    plt.scatter(stellar_mass, manga_star, s=1, color='navy', label='rmins')
    plt.scatter(stellar_mass, mpa_sfr_manga, s=1, color='lightsteelblue', label='MPA-JHU')
    plt.scatter(stellar_mass, manga_starb, s=1, color='cyan', label='MPA-JHU')
    plt.xlim(0)
    plt.ylim(-1.6, 0.1)
    plt.xlabel('Stellar Mass (M☉/h^2)')
    plt.ylabel('Log(Star Formation Rate)')
    plt.show()
    pdf.savefig()




