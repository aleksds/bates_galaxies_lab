# Aleks Diamond-Stanic 20170706
# based on code written by John Moustakas

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#import fitsio
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

massdir = os.getenv('ISEDFITDIR')
photfile = 'sg_fluxtable_nm.txt'
isedfile = os.path.join(massdir, 'massprofiles_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits.gz')
kcorrfile = os.path.join(massdir, 'massprofiles_fsps_v2.4_miles_chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz')

print('Reading {}'.format(photfile))
phot = ascii.read(photfile)
phot[:2]

