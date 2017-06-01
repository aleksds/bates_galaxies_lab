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

# define a function to plot an absorption line profile
def plot_profile(wave, flux, line, name, fosc):
    # define the velocity scale [AA / s]
    vel_aas = (wave - line*(1+zem[i])) / (line*(1+zem[i])) * c
    # convert to [km / s]
    vel_kms = vel_aas.to('km/s')
    # define parameters for the x and y axes
    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(axis='x', labelsize=xls)
    ymax = 3.e12
    ax.set_ylim(0., ymax)
    # make the plot (using equations 5 and 8 from Savage & Sembach 1991)
    ax.plot(vel_kms, np.log(1/flux) / 2.654e-15 / (fosc * line))
    # include the name of the line
    plt.text(xmin+0.03*(xmax-xmin), ymax*0.6, name)
    # mark the approximate centroid velocity
    plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
    plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    # label other lines for clarity
    for k in range(0, len(lines)):
        vel_off_aas = (lines[k] - line) / line * c
        vel_off_kms = vel_off_aas.to('km/s') / (u.km / u.s)
        plt.axvline(x=vcen[i]+vel_off_kms, ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
            
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

# set relevant parameters for the figure
filename = 'all_nh.pdf'
xls = 7.
minorLocator = AutoMinorLocator()

with PdfPages(filename) as pdf:

    for i in range(0, len(gal)):

        # read in the spectrum
        datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(500. * u.km * u.s)

        fig = plt.figure()
        plt.suptitle(gal[i])

        # loop over each spectral line
        for j in range(0, len(lines)):
            ax = fig.add_subplot(4,2,j+1)
            plot_profile(wave, flux, lines[j], names[j], fosc[j])
        
        pdf.savefig()
        plt.close()
    
    os.system("open %s &" % filename)
