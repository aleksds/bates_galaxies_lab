# Rebecca Minsley
# April 29th 2019

import os
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
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


mpl8_dir = os.environ['MANGADIR_MPL8']  # directory should be saved in .bash_profile

drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')            # read in information from drpall file
nsa = fits.open(mpl8_dir + '1-nsa_v1_0_1.fits')             # read in information from NSA catalog
gal = fits.open(mpl8_dir + 'gal_info_dr7_v5_2.fit')         # MPA-JHU catalog galaxy data
sfr = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')      # MPA-JHU catalog star-formation
mass = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')     # MPA-JHU catalog galactic mass

drpdata = drp[1].data       # drpall
nsa_data = nsa[1].data      # nsa catalog
galdata = gal[1].data       # MPA-Jhu
sfrdata = sfr[1].data       # MPA-Jhu
massdata = mass[1].data     # MPA-Jhu


m_survey = drpdata.field('srvymode') == 'MaNGA dither'      # boolian -> True at MaNGA Data
stellar_mass = drpdata.field('nsa_sersic_mass')[m_survey]   # MaNGA gal stellar mass
ba = drpdata.field('nsa_sersic_ba')[m_survey]               # sersic profile of Manga Galaxies
plateifu = drpdata.field('plateifu')[m_survey]              # MaNGA plateifu
manga_ra = drpdata.field('objra')[m_survey]                 # MaNGA declination
manga_dec = drpdata.field('objdec')[m_survey]

nsa_ra = nsa_data.field('RA')       # nsa ra
nsa_dec = nsa_data.field('DEC')     # nsa dec



# match the two catalogs on the basis of right ascension and declination
c = SkyCoord(ra=manga_ra*u.degree, dec=manga_dec*u.degree)
catalog = SkyCoord(ra=nsa_ra*u.degree, dec=nsa_dec*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)
# match the MPA-JHU and MaNGA catalogs on the basis of RA and DEC
mpa_ra = galdata.field('RA')
mpa_dec = galdata.field('DEC')
bad = mpa_dec < -90
mpa_dec[bad] = -90
mpa_cat = SkyCoord(ra=mpa_ra*u.degree, dec=mpa_dec*u.degree)
idx_mpa, d2d_mpa, d3d_mpa = c.match_to_catalog_sky(mpa_cat)



# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
