# Aleks Diamond-Stanic
# 20161004
#
# The main goal of this code is to make plots of the apparent optical
# depth for relevant Mg and Fe absorption lines.  Based on hires_plot_abs.py
#
# There is one galaxy per page, and there are 4 rows x 2 columns = 8
# panels with the following lines: Mg I 2852, Mg II 2803, Mg II 2796,
# Fe II 2600, Fe II 2586, Fe II 2382, Fe II 2374, Fe II 2344
#

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator

# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')

# wavelengths of relevant absorption lines
mgi2852 = 2852.96328 * u.AA
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
feii2600 = 2600.1724835 * u.AA
feii2586 = 2586.6495659 * u.AA
feii2382 = 2382.7641781 * u.AA
feii2374 = 2374.4603294 * u.AA
feii2344 = 2344.2129601 * u.AA
 
#oscillator strengths (f_ij)
f2852 = 1.83
f2803 = 0.3058
f2796 = 0.6155
f2600 = 0.2394
f2586 = 0.069126
f2382 = 0.320
f2374 = 0.0313
f2344 = 0.1142


# arrays of line wavelengths and names
lines = [mgii2796, mgii2803, feii2586, feii2600, feii2374, feii2382, feii2344, mgi2852]
fosc = [f2796, f2803, f2586, f2600, f2374, f2382, f2344, f2852]
names = ['Mg II 2796', 'Mg II 2803', 'Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344', 'Mg I 2852']

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']

# Calculate ratio of MgI, MgII to FeII

#mass of electron
mass_e = 9.11**(-31)

# tau is related to 1/e and inversely proportional to optical depth
def calculate(wave, flux, line, fosc):
    # define the velocity scale [AA / s]
    vel_aas = (wave - line*(1+zem[i])) / (line*(1+zem[i])) * c
    # convert to [km / s]
    vel_kms = vel_aas.to('km/s')
    # label other lines for clarity

    tau = np.log(1/flux) 

    for k in range(0, len(lines)):
        vel_off_aas = (lines[k] - line) / line * c
        vel_off_kms = vel_off_aas.to('km/s') / (u.km / u.s)
    return(tau)  

for i in range(0, len(gal)):
    # read in the spectrum
    datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
    data = ascii.read(datafile)
    wave = data['wv'] * u.AA
    flux = data['norm']  #intensity of the light you obseved   ;  normalized means you check observed light from expected light

success = np.array([])
# loop over each spectral line
for j in range(0, len(lines)):
    blah = calculate(wave, flux, lines[j], fosc[j])
    success = np.append(success,blah)

#wavelength to be squared!
wavelength = np.array([2852.96328, 2803.5314853, 2796.3542699, 2600.1724835, 2586.6495659, 2382.7641781, 2374.4603294, 2344.2129601])  

#absorption oscillator strength
f_ij = np.array ([1.83, 0.3058, 0.6155, 0.2394, .069126, 0.320, 0.0313, 0.1142])

for i in range(0, len(f_ij)):
    col_dens = mass_e * c * success[i] / (np.pi * np.exp(2) * f_ij[i] * wavelength**2)
    print(names[i],col_dens[i])



