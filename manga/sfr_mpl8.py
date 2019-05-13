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


# def drp_name(plate_ifu):
#     for i in range(len(plate_ifu)):
#         p, fu = good_plates[i].split('-')
#         file_name = mpl8_dir + 'HYB10-MILESHC-MILESHC/' + str(p) + '/' + str(fu) + '/manga-' + str(p) + '-' + str(fu) \
#                     + '-MAPS-HYB10-MILESHC-MILESHC.fits.gz'
#
#     return file_name



cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)      # define a standard cosmology


mpl8_dir = os.environ['MANGADIR_MPL8']  # directory should be saved in .bash_profile

drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')            # read in information from drpall file
dap = fits.open(mpl8_dir + 'dapall-v2_5_3-2.3.0.fits')      # read in information from dapall file
nsa = fits.open(mpl8_dir + '1-nsa_v1_0_1.fits')             # read in information from NSA catalog

gal = fits.open(mpl8_dir + 'gal_info_dr7_v5_2.fit')         # MPA-JHU catalog galaxy data
sfr = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')      # MPA-JHU catalog star-formation
mass = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')     # MPA-JHU catalog galactic mass


drpdata = drp[1].data       # drpall
dapdata = dap[1].data       # dapall
nsa_data = nsa[1].data      # nsa catalog

galdata = gal[1].data       # MPA-Jhu
sfrdata = sfr[1].data       # MPA-Jhu
massdata = mass[1].data     # MPA-Jhu

# drp
m_survey = drpdata.field('srvymode') == 'MaNGA dither'      # boolean -> True at MaNGA Data
stellar_mass = drpdata.field('nsa_sersic_mass')             # MaNGA gal stellar mass
ba = drpdata.field('nsa_sersic_ba')[m_survey]               # sersic profiles
plateifu = drpdata.field('plateifu')[m_survey]              # MaNGA plateifu
manga_ra = drpdata.field('objra')                           # declination
manga_dec = drpdata.field('objdec')

# dap
sfr = dapdata.field('SFR_TOT')                    # sfr rate h_flux lumin from DAP file


# nsa catalog
nsa_ra = nsa_data.field('RA')               # nsa ra
nsa_dec = nsa_data.field('DEC')             # nsa dec
th50 = nsa_data.field('ELPETRO_TH50_R')     # th50 concentration array
th90 = nsa_data.field('ELPETRO_TH90_R')     # th90 concentration array


# match the two catalogs on the basis of right ascension and declination
c = SkyCoord(ra=manga_ra*u.degree, dec=manga_dec*u.degree)
catalog = SkyCoord(ra=nsa_ra*u.degree, dec=nsa_dec*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)
th50_manga = th50[idx]
th90_manga = th90[idx]
c_manga = (th90_manga/th50_manga)[m_survey]                 # MaNGA galaxy concentrations
late = c_manga < 2.6                                        # late-type concentrations (<2.6)
edge = ba[late] < 0.3                                       # Defining edge-on Galaxies
dist = drpdata.field('nsa_zdist')[m_survey][late][edge]     # multiply by c/H0 for Mpc.

# match the MPA-JHU and MaNGA catalogs on the basis of RA and DEC
mpa_ra = galdata.field('RA')
mpa_dec = galdata.field('DEC')
bad = mpa_dec < -90
mpa_dec[bad] = -90
mpa_cat = SkyCoord(ra=mpa_ra*u.degree, dec=mpa_dec*u.degree)
idx_mpa, d2d_mpa, d3d_mpa = c.match_to_catalog_sky(mpa_cat)


# MPA-JHU catalog matched data
mpa_sfr = sfrdata.field('MEDIAN')[idx_mpa]
mpa_mass = massdata.field('MEDIAN')[idx_mpa]
mpa_sfr_manga = mpa_sfr[m_survey][late][edge]       # MPA-JHU late-type edge-on MaNGA galaxy sfr
mpa_mass_manga = mpa_mass[m_survey][late][edge]     # MPA-JHU late-type edge-on MaNGA galaxy stellar mass


good_plates = plateifu[late][edge]       # late-type edge-on MaNGA plateifus
manga_star = np.zeros(len(good_plates))
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

        z = (dist[i] * (const.c / cosmo.H0.to(u.m/(u.m * u.s)))).to(u.cm)
        obs_ha_flux = (ha_flux * (10 ** (-0.4 * k_halpha * color_excess)))    # in 10^-17 erg/s/cm2/spaxel
        obs_ha_flux = obs_ha_flux * 10**-17 * u.erg/u.s/(u.cm ** 2)

        lum = 4 * math.pi * z**2 * obs_ha_flux
        log_sfr = math.log(lum.value, 10) - 41.27

        manga_star[i] = log_sfr

        # mpl8_sfr = dap_data['SFR_TOT']
