#This program loops through all the galaxies we are focused on
#2796 and 2803 Mg 
#2600 and 2382 Fe
#Creates a page with four graphs for each galaxy

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d


# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
feii2600 = 2600.1724835 * u.AA
feii2382 = 2382.7641781 * u.AA
 
#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
f2600 = 0.2394
f2382 = 0.320

# array of names, lines of interest, oscillator strength:
names = ['Mg II 2796', 'Fe II 2600']
lines = [mgii2796, feii2600]
fosc = [f2796, f2600]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']

def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
#if len(cor_col) < 1352:
        #cor_col = np.append(cor_col, cor_col[1350])
    return cor_col

def velocity(vel, col_dens):
    cor_vel = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_vel = np.append(cor_vel, vel[i])
    #if len(cor_vel) < 1352:
     #   cor_vel = np.append(cor_vel, cor_vel[1350])
    return cor_vel

Mg = 0
Fe = 1

counter1 = 0
minorLocator = AutoMinorLocator()
filename = 'Ratio.pdf'
with PdfPages(filename) as pdf:
    datafile = dir+gal[2]+'/'+gal[2]+'_stitched_v1.txt'
    data = ascii.read(datafile)
    wave = data['wv'] * u.AA
    flux = data['norm'] 
    vel_kms = np.zeros([len(lines),len(wave)])
    # define the velocity scale [km / s]
    for i in range(0, len(lines)):
        vel_kms[i] = ((wave-lines[i]*(1+zem[2]))/(lines[i]*(1+zem[2]))) * c_kms
    tau = np.zeros([len(lines),len(flux)])
    # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
    for j in range(0, len(lines)):
        blah = np.log(1/flux)
        tau[j] = blah
        
    col_dens = np.zeros([len(lines), len(flux)])
    #calculating column density
    for i in range(0, len(fosc)):
        col_dens[i] = tau[i] / (2.654E-15 * fosc[i] * (wave/(1+zem[2])) * fosc[i])

    #Mg is 0
    #Fe is 1
    fig = plt.figure()
    cor_col_Fe = column(vel_kms[1], col_dens[1])
    cor_col_Mg = column(vel_kms[0], col_dens[0])
    cor_vel_Fe = velocity(vel_kms[1], col_dens[1])
    cor_vel_Mg = velocity(vel_kms[0], col_dens[0])

    
    
    #INTERPOLATE:
    f0 = interp1d(vel_kms[0], flux)
    f1 = interp1d(vel_kms[1], flux)

    vel_new0 = np.linspace(-3000, 500, num=3500, endpoint=True)
    #vel_new1 = np.linspace(-3000, 500, num=len(flux), endpoint=True)

    flux_king = f0(vel_new0)
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(vel_kms[0], flux, marker='+')
    ax.plot(vel_new0, flux_king)
    plt.ylabel('Flux')
    ax.set_xlim(-3000,-1500)
    ax.set_ylim(0, 2)
    
    pdf.savefig()
    plt.close()
            
os.system("open %s &" % filename)
