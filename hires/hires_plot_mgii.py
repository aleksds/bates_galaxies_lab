# Aleks Diamond-Stanic
# 20160517 --> 20160926
#
# the goals of this code include the following:
#
# (1) make plots that visualize the spectra in the wavelength region
#     surrounding Mg II
#
# (2) show the line profiles for the 2796 and 2803 transitions
#     separately for each of 14 galaxies
#
# (3) visualize the velocity structure of the gas traced by MgII using
#     the 2796 lines for the highest velocities and the 2803 for
#     velocities closer to v=0

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator

# define a function to set plot parameters
def plot_setup():
    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.set_ylim(0., 1.5)
    ax.tick_params(axis='x', labelsize=xls)
    ax.tick_params(axis='y', labelsize=yls)

# define a function to mark the centroid velocity
def plot_vline():
    plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
    plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')

# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA

# arrays of galaxy names, redshifts, approximate centroid velocities,
# approximate maximum velocities, and velocities to "flip" between the 2796
# and 2803 line when constructing an empirical mgii velocity profile
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip'] * u.km / u.s

# sort by maximum outflow velocity
vmax_index = np.flipud(np.argsort(vmax))

# define velocity ranges to plot the profile
vb = -3500. * u.km / u.s
vr = 1500. * u.km / u.s

# make plot showing one galaxy per page
filename = 'all_mgii.pdf'
xls = 7.
yls = 10.
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

        # define the 2803 and 2796 velocity scales
        vel_mgii_2803_aa = (wave - mgii2803*(1+zem[i])) / (mgii2803*(1+zem[i])) * c
        vel_mgii_2803 = vel_mgii_2803_aa.to('km/s')
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[i])) / (mgii2796*(1+zem[i])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # show the profile in 2796 velocity units
        ax = fig.add_subplot(4,1,1)
        plot_setup()
        ax.plot(vel_mgii_2796, flux)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796', color='blue')
        plot_vline()

        # show the profile in 2803 velocity units
        ax = fig.add_subplot(4,1,2)
        plot_setup()
        ax.plot(vel_mgii_2803, flux, color='red')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2803', color='black')
        plot_vline()

        # show the profile in both 2796 and 2803 velocity units
        ax = fig.add_subplot(4,1,3)
        plot_setup()
        ax.plot(vel_mgii_2803, flux, color='red')
        ax.plot(vel_mgii_2796, flux, color='blue')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796+2803')
        plot_vline()

        # show the velocity profile using the 2796 line on the blue
        # side and the 2803 line on the red side
        g2796 = (vel_mgii_2796 > vb) & (vel_mgii_2796 < vflip[i])
        g2803 = (vel_mgii_2803 > vflip[i]) & (vel_mgii_2796 < vr)
        ax = fig.add_subplot(4,1,4)
        plot_setup()
        ax.plot(vel_mgii_2803[g2803], flux[g2803], color='red')
        ax.plot(vel_mgii_2796[g2796], flux[g2796], color='blue')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796+2803')
        plot_vline()
        
        pdf.savefig()
        plt.close()
    
    os.system("open %s &" % filename)


# make plot showing all galaxies on each page
filename = 'all_mgii_2page.pdf'
xls = 8.
yls = 8.


with PdfPages(filename) as pdf:

    fig = plt.figure()
    
    for i in range(0, len(gal)):

        # read in the spectrum, sorted by outflow velocity
        indx = vmax_index[i]
        datafile = dir+gal[indx]+'/'+gal[indx]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(1500. * u.km * u.s)

        # define the 2796 velocity scale
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[indx])) / (mgii2796*(1+zem[indx])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # plot the profiles using the 2796 velocity scale
        ax = fig.add_subplot(5,3,i+1)
        plot_setup()
        ax.plot(vel_mgii_2796, flux, color='black', linewidth=0.5)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, gal[indx])

    pdf.savefig()
    plt.close()
    
    fig = plt.figure()
        
    for i in range(0, len(gal)):

        # read in the spectrum, sorted by outflow velocity
        indx = vmax_index[i]
        datafile = dir+gal[indx]+'/'+gal[indx]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(1500. * u.km * u.s)

        # define the 2803 and 2796 velocity scale
        vel_mgii_2803_aa = (wave - mgii2803*(1+zem[indx])) / (mgii2803*(1+zem[indx])) * c
        vel_mgii_2803 = vel_mgii_2803_aa.to('km/s')
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[indx])) / (mgii2796*(1+zem[indx])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_mgii_2796 > vb) & (vel_mgii_2796 < vflip[indx])
        g2803 = (vel_mgii_2803 > vflip[indx]) & (vel_mgii_2803 < vr)

        # plot the profiles using the 2796 profile on the blue side
        # and the 2803 profile on the red side
        ax = fig.add_subplot(5,3,i+1)
        plot_setup()
        ax.plot(vel_mgii_2803[g2803], flux[g2803], color='red', linewidth=0.5)
        ax.plot(vel_mgii_2796[g2796], flux[g2796], color='blue', linewidth=0.5)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, gal[indx])
        
    pdf.savefig()
    plt.close()
    
    os.system("open %s &" % filename)
