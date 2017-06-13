# Aleks Diamond-Stanic, 20170606
# code to read in and plot an SDSS spectrum

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy import constants as const
from astropy import units as u

dir = os.environ['SDSSDIR']
file = dir + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)
hdulist.info()
coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']

flux = hdulist[1].data


npix = len(flux)
index = np.arange(npix)
wavelength = 10.**(coeff0 + coeff1*index)

# speed of light in Angstroms per second
c = c = const.c.to('AA/s')
minorLocator = AutoMinorLocator()

wavelength_of_interest = ('Fe2600', 'MgII2796', 'MgII2803', 'MgI2852')
z = .60329
rest_wavelength = np.array([2600, 2796, 2803, 2852])




# Observed Wavelength
obs_wavelength = np.array([])
for i in range(0,len(rest_wavelength)):
    wavelength_calc = rest_wavelength[i] * (1 + z)
    obs_wavelength = np.append(obs_wavelength, wavelength_calc)

#Observed Velocity
for i in range(0,len(rest_wavelength)):
    velocity_calc = (obs_wavelength- rest_wavelength[i]) * c / rest_wavelength[i]
obs_velocity = np.array([velocity_calc.to('km/s')])

#velocity conversion
for i in range(0, len(wavelength)):    
    vel_big = ((wavelength*(1+z)) - wavelength[i]) * c / wavelength[i]

# convert to [km / s]
vel_big_graph = np.array([vel_big.to('km/s')])
    
# define a function to plot an absorption line profile
def sdss_plot_wavelength(wavelengths, name):

    # define parameters for the x and y axes
    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(axis='x', labelsize=7)
    ax.set_ylim(0., 15.)
    # make the plot
    ax.plot(wavelength, flux)
    # include the name of the line
    plt.text(xmin+0.03*(xmax-xmin), 0.15, name)
    # mark the approximate centroid velocity
    plt.axvline(x=wavelengths, ymin=0., ymax = 15., linewidth=1, color='k', linestyle='dotted')
    plt.axvline(x=wavelengths + 30., ymin=0., ymax = 15., linewidth=0.5, color='k')
    plt.axvline(x=wavelengths - 30., ymin=0., ymax = 15., linewidth=0.5, color='k')


def sdss_plot_velocity(wave, interest_wave, flux, name, target):

    # define the velocity scale
    vel_aas = (wave/(1+z) - interest_wave) * c / interest_wave 
    # convert to [km / s]
    vel_kms = vel_aas.to('km/s')
    #define parameters for the x and y axes
    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(axis='x', labelsize=7)
    ax.set_ylim(0., 15.)
    # make the plot
    ax.plot(vel_kms, flux)
    # include the name of the line
    plt.text(xmin+0.03*(xmax-xmin), 0.15, name)
    # mark the approximate velocities of interest
    plt.axvline(x=target, ymin=0., ymax = 15., linewidth=1, color='k', linestyle='dotted')
    plt.axvline(x=target +30., ymin=0., ymax = 15., linewidth=0.5, color='k')
    plt.axvline(x=target -30., ymin=0., ymax = 15., linewidth=0.5, color='k')


filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:


#Big Graph


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(wavelength, flux)
    ax.set_ylim(0., 15.)
    for i in range(0, len(obs_wavelength)):
        # mark the observed wavelengths of interest
        plt.axvline(x=obs_wavelength[i], ymin=0., ymax = 15., linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=obs_wavelength[i]+30., ymin=0., ymax = 15., linewidth=0.5, color='k')
        plt.axvline(x=obs_wavelength[i]-30., ymin=0., ymax = 15., linewidth=0.5, color='k')
    plt.xlabel('Wavelength',fontsize=14)
    plt.ylabel('Flux',fontsize=14)
    pdf.savefig()
    plt.close()


#Big Graph Based on Velocity

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(vel_big_graph[0], flux)
    ax.set_ylim(0., 15.)
    for i in range(0, len(obs_velocity)):
        # mark the observed wavelengths of interest
        plt.axvline(x=obs_velocity[i].all(), ymin=0., ymax = 15., linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=(obs_velocity[i]+30.).all(), ymin=0., ymax = 15., linewidth=0.5, color='k')
        plt.axvline(x=(obs_velocity[i]-30.).all(), ymin=0., ymax = 15., linewidth=0.5, color='k')
    plt.xlabel('Velocity (observed)',fontsize=14)
    plt.ylabel('Flux',fontsize=14)
    pdf.savefig()
    plt.close()
    
    
#Graphs of Wavelengths of Interest

    for i in range(0, len(wavelength_of_interest)):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        xmin = obs_wavelength[i] - 100
        xmax = obs_wavelength[i] + 100
        ax.plot(wavelength, flux)
        plt.xlabel('Wavelength',fontsize=14)
        plt.ylabel('Flux',fontsize=14)
        ax.set_ylim(0., 15.)
        plt.suptitle(wavelength_of_interest[i])
        sdss_plot_wavelength(obs_wavelength[i], wavelength_of_interest[i])

        pdf.savefig()
        plt.close()

#Graphs of the Velocities of the Wavelengths of Interest

    for i in range(0, len(wavelength_of_interest)):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        xmin = -3000
        xmax = 500
        
        plt.xlabel('Velocity km/s (observed)',fontsize=14)
        plt.ylabel('Flux',fontsize=14)

        plt.suptitle(wavelength_of_interest[i])
        sdss_plot_velocity(wavelength,rest_wavelength[i], flux, wavelength_of_interest[i], -1600.)

        pdf.savefig()
        plt.close()

    os.system("open %s &" % filename)

