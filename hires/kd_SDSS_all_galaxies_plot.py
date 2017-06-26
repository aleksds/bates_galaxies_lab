#Kwamo, 20170615
#Read all SDSS spectrums and make a plot of flux versus wavelength, velocity
#
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
import glob


# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

#data directory for SDSS data
dir = os.environ['SDSSDIR']
file = glob.glob(dir+'*.fits')
#read in the relevant sdss files
for i in range(0, len(file)):
    hdulist = fits.open(file[i])
flux = hdulist[1].data
coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']
npix = len(flux)
index = np.arange(npix)
wave = np.array([10.**(coeff0 + coeff1*index)])
## SDSS PLOT WORK!!

#redshift coefficient
z = 0.60329
#wavelengths:
w_int = ['MgII2796','MgII2803', 'MgI2852', 'FeII2600']
w_rest = np.array([2796, 2803 , 2852, 2600])
w_obs = []
#define w_obs:
for i in range(0, len(w_rest)):
    temp = w_rest[i]*(1+z)
    w_obs.append(temp)
print(w_obs)

#wavelength to velocity conversion:
vel_obs = np.array([])
vel_aas = np.zeros([len(w_obs),len(flux)])
vel_aas_0 = []
vel_aas_1 = []
vel_aas_2 = []
vel_aas_3 = []

vel_aas = [[vel_aas_0], [vel_aas_1], [vel_aas_2], [vel_aas_3]]
vel_kms = np.zeros([len(w_obs),len(flux)])

# define the velocity scale [AA / s]

for i in range(0, len(w_obs)):

    vel_aas[i] = ((wave[0]-w_obs[i])/w_obs[i]) * c
    # convert to [km / s]
    vel_kms[i] = np.array([vel_aas[i].to('km/s')])
    vel_kms[i] = vel_aas[i]

for i in range(0, len(w_obs)):
    vel_temp = (w_obs[i]) / (w_obs[i]*(1+z)) * c
    vel_kms_temp = np.array([vel_temp.to('km/s')], dtype = object)
    vel_obs = np.append(vel_obs, vel_kms_temp)
    
# #Question: how long should tau be?  1 value for each line, same # as flux for each line, just the same # as flux in total?
# tau = np.zeros([len(w_obs),len(flux)])
# # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
# for j in range(0, len(w_obs)):
#     boom = np.log(1/flux)
#     tau[j] = boom

    #  NEED TO FIND THE SDSS VERSION OF FOSC (OSCILLATOR STRENGTH) TO FINISH SDSS COL_DENS CALCULATION
#Column Density Calc
#col_dens[i] = tau[i] / (2.654E-15 * fosc[i] * (wave/(1+z)) * fosc[i])


# SDSS PLOTTING

minorLocator = AutoMinorLocator()
filename = 'J0826_sdss.pdf'

with PdfPages(filename) as pdf:

#prints the first graph 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(wave[0], flux)

    for i in range(0, len(w_rest)):

        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    ax.set_ylim(0., 15.)
    ax.set_xlim(3800,9000)
    pdf.savefig()
    plt.close()

    #prints each of the zoomed in graphs

    for i in range(0, len(w_int)):
        # read in the spectrum
        fig2 = plt.figure()
        ax = fig2.add_subplot(1,1,1)
        ax.plot([(1, 2), (3, 4)], [(4, 3), (2, 3)])
        plt.title(w_int[i])
        ax.set_xlim(w_obs[i]-200, w_obs[i]+200)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.plot(wave[0], flux)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+5., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-5., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        pdf.savefig()
        plt.close()
        
    # #prints the graph with velocity
    # for i in range(0, len(w_int)):
    #     fig3 = plt.figure()
    #     ax = fig3.add_subplot(1,1,1)
    #     ax.plot(vel_kms[i], flux)
    #     plt.title(w_int[i])
    #     plt.axvline(x=-1600, ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
    #     plt.axvline(x=-1600+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    #     plt.axvline(x=-1600-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    #     plt.xlabel('Velocity (Observed)')
    #     plt.ylabel('Flux')
    #     ax.set_ylim(0., 15.)
    #     ax.set_xlim(-3000,500)
    #     pdf.savefig()
    #     plt.close()
            
            
    os.system("open %s &" % filename)
