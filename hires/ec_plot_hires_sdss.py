#Eve Cinquino, 20170606
#Read in an SDSS spectrum and make a plot of flux versus wavelength
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

#data directory for SDSS data
dir = os.environ['SDSSDIR']
file = dir + 'J1713.fits'
hdulist = fits.open(file)

# define values for SDSS data
coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']
flux = hdulist[1].data
model = hdulist[1].data['model']
ivar = hdulist[1].data['ivar']
npix = len(flux)
index = np.arange(npix)
wavelength = np.array([10.**(coeff0 + coeff1*index)])

flux1 = np.zeros([len(flux)])
for i in range(0, len(flux)):
    flux1[i] = flux[i][0]/model[i]

dir2 = os.environ['HIRESDIR']
# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir2+'gal_info.txt')
gal = gal_info['gal']
vcen = gal_info['vcen']
datafile = dir2+gal[0]+'/'+gal[0]+'_stitched_v1.txt'
print(datafile)
data = ascii.read(datafile)
wave2 = np.array([data['wv'] * u.AA])
flux2 = data['norm']
# set relevant parameters for the figure
filename2 = 'all_abs.pdf'



# speed of light in Angstroms per second
c = const.c.to('km/s')
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
vel_aas = np.zeros([len(w_obs),3842])

#vel_aas = [[vel_aas_0], [vel_aas_1], [vel_aas_2], [vel_aas_3]]
vel_kms = np.zeros([len(w_obs),3842])
# define the velocity scale [AA / s]
for i in range(0, len(w_obs)):
    vel_aas[i] = ((wavelength[0]-w_obs[i])/w_obs[i]) * c
    # convert to [km / s]
    #vel_kms[i] = np.array([vel_aas[i].to('km/s')])
    vel_kms[i] = vel_aas[i]

vel_kms2 = np.zeros([len(w_obs), len(flux2)])
for i in range(0, len(w_obs)):
    vel_kms2[i] = ((wave2[0]-w_obs[i])/w_obs[i]) * c

for i in range(0, len(w_obs)):
    vel_temp = (w_obs[i]) / (w_obs[i]*(1+z)) * c
    vel_kms_temp = np.array([vel_temp.to('km/s')], dtype = object)
    vel_obs = np.append(vel_obs, vel_kms_temp)

#PLOTTING
minorLocator = AutoMinorLocator()
filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:
    #prints the first graph 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #ax.plot(wavelength[0], flux)
    plt.plot(wave2[0], flux2, color = 'r', label='HIRES', linewidth=0.5)
    plt.plot(wavelength[0], flux1, color = 'g', label='SDSS', linewidth=1)
    
    #plt.legend(loc='best')
    #plt.show()
    for i in range(0, len(w_rest)):
        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-30., ymin=0., ymax = 1.5, linewidth=1, color='k')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    ax.set_ylim(0., 2.)
    ax.set_xlim(3800,9000)
    pdf.savefig()
    plt.close()

    #prints each of the zoomed in graphs
    for i in range(0, len(w_int)):
        # read in the spectrum
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        #ax.plot([(1, 2), (3, 4)], [(4, 3), (2, 3)])
        plt.title(w_int[i])
        ax.set_xlim(w_obs[i]-100, w_obs[i]+100)
        ax.set_ylim(-1, 2)
        #ax.xaxis.set_minor_locator(minorLocator)
        plt.plot(wave2[0], flux2, color = 'r', label='HIRES', linewidth=0.5)
        plt.plot(wavelength[0], flux1, color = 'g', label = 'SDSS', linewidth=1)
        
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.axvline(x=w_obs[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs[i]+5., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs[i]-5., ymin=0., ymax = 1.5, linewidth=1, color='k')

    #prints the graph with velocity
        ax = fig.add_subplot(3,1,3)
        plt.plot(vel_kms2[i], flux2, color = 'r', label='HIRES', linewidth=0.5)
        plt.plot(vel_kms[i], flux1, linewidth=1)
        plt.axvline(x=-1600, ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=-1600+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=-1600-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.xlabel('Velocity (Observed)')
        plt.ylabel('Flux')
        ax.set_ylim(0., 2.)
        ax.set_xlim(-3000,500)
        pdf.savefig()
        plt.close()
    
    os.system("open %s &" % filename)
