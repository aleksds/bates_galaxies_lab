# Jose Ruiz

# Goal for the code:
# Make an unbinned plot of Cover Fraction vs velocity next to a
#plot of Flux vs Velocity, for the a velocity range of -1500 to -1000 km/s


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
vb = -100. * u.km / u.s
vr = 0. * u.km / u.s


# Make plot showing velocity profile and Covering Fraction Profile
filename='Test02_mgii.pdf'
xls=5.
yls=5.
minorLocator = AutoMinorLocator()
with PdfPages(filename) as pdf:
    for i in range(0,len(gal)):

        #read in the spectrum
        datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-1500. * u.km * u.s)
        xmax = np.array(-1000. * u.km * u.s)
        Intensity=1-flux

        fig = plt.figure()
        plt.suptitle(gal[i])

         # define the 2803 and 2796 velocity scales
        vel_mgii_2803_aa = (wave - mgii2803*(1+zem[i])) / (mgii2803*(1+zem[i])) * c
        vel_mgii_2803 = vel_mgii_2803_aa.to('km/s')
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[i])) / (mgii2796*(1+zem[i])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

         # show the profile in 2796 velocity units, for Vel profile
        ax = fig.add_subplot(4,1,1)
        plot_setup()
        ax.set_ylim(0., 1.5)
        ax.plot(vel_mgii_2796, flux)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796', color='blue')
        plot_vline()
        plt.ylabel("Flux")

        # show the profile in 2803 velocity units, for Vel profile
        ax = fig.add_subplot(4,1,2)
        plot_setup()
        ax.set_ylim(0., 1.5)
        ax.plot(vel_mgii_2803, flux, color='red')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2803', color='black')
        plot_vline()
        plt.ylabel("Flux")

         # show the profile in 2796 velocity units, for Cf profile
        ax = fig.add_subplot(4,1,3)
        plot_setup()
        ax.set_ylim(0., 1.0)
        ax.plot(vel_mgii_2796, Intensity)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796', color='blue')
        plot_vline()
        plt.ylabel("$C_f$")

        # show the profile in 2803 velocity units, for Cf profile
        ax = fig.add_subplot(4,1,4)
        plot_setup()
        ax.set_ylim(0., 1.0)
        ax.plot(vel_mgii_2803, Intensity, color='red')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2803', color='black')
        plot_vline()
        plt.ylabel("$C_f$")
        plt.xlabel("Velocity [km/s]")

        pdf.savefig()
        plt.close()

    os.system("open %s &" % filename)

      
        
