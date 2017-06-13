#Eve Cinquino
#20170607
#
#The goal of this project is to determine the ratio of Mg to Fe and to then make plots of their the apparent optical depth for relevant absorption lines
#
#First: define important variables, wavelengths, names, c, fosc, gal_info
#Second: define functions, need to find tau, column density, then ratio
#Print out the ratio
#
#Then: make graphs
#
#SOMETHING I IGNORED: THERE ARE 14 DIFFERENT GALAXIES, I THINK MY CODE ONLY RUNS THROUGH J0826


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
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

# wavelengths of relevant absorption lines
mgi2852 = 2852.96328 * u.AA
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
feii2600 = 2600.1724835 * u.AA
feii2586 = 2586.6495659 * u.AA
feii2382 = 2382.7641781 * u.AA
feii2374 = 2374.4603294 * u.AA
feii2344 = 2344.2129601 * u.AA
 
#oscillator strengths (fosc)
f2852 = 1.83
f2803 = 0.3058
f2796 = 0.6155
f2600 = 0.2394
f2586 = 0.069126
f2382 = 0.320
f2374 = 0.0313
f2344 = 0.1142

# array of names, lines of interest, oscillator strength:
names = ['Mg II 2796', 'Mg II 2803', 'Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344', 'Mg I 2852']
lines = [mgii2796, mgii2803, feii2586, feii2600, feii2374, feii2382, feii2344, mgi2852]
fosc = [f2796, f2803, f2586, f2600, f2374, f2382, f2344, f2852]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']

for i in range(0, len(gal)):
    # read in the spectrum
    datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
    data = ascii.read(datafile)
    wave = data['wv'] * u.AA
    flux = data['norm']  #intensity of the light you obseved   ;  normalized means you check observed light from expected light

vel_kms = np.zeros([len(lines),50515])
# define the velocity scale [km / s]
for i in range(0, len(lines)):
    vel_kms[i] = ((wave-lines[i]*(1+zem[0]))/(lines[i]*(1+zem[0]))) * c_kms

#Question: how long should tau be?  1 value for each line, same # as flux for each line, just the same # as flux in total?
tau = np.zeros([len(lines),len(flux)])
# loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
for j in range(0, len(lines)):
    blah = np.log(1/flux)
    tau[j] = blah

col_dens = np.zeros([len(lines), len(flux)])
#calculating column dens, again how many values should this be?? 
for i in range(0, len(fosc)):
    #col_dens[i] = mass_e * c * tau[i] / (np.pi * np.exp(2) * fosc[i] * wave**2)
    col_dens[i] = tau[i] / (2.654E-15 * fosc[i] * (wave/(1+zem[0])) * fosc[i]) 
#now to find the ratio???  There are 15 ratios to find, that compare Mg to Fe
#index values for above arrays for Mg and Fe:
Mg = [0, 1, 7]
Fe = [2, 3, 4, 5, 6]
Ratio = np.zeros([15, len(col_dens[0])])
counter = 0
for i in range(0, len(Mg)):
    for j in range(0, len(Fe)):
         Ratio[counter] = col_dens[Mg[i]]/col_dens[Fe[j]]
         counter += 1
#Successfully made a 15 by 50515 array for ratio, comparing each Mg line to each Fe line

#Now, GRAPHING: (vel/wave vs. ratio)
minorLocator = AutoMinorLocator()
filename = 'Ratio.pdf'
with PdfPages(filename) as pdf:
    for i in range(0, len(Mg)):
        for j in range(0, len(Fe)):
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.set_xlim(-3000,500)
            ymax = 3.e12
            ax.set_ylim(0., ymax)
            plt.title("Comparing %s and %s" % (names[Mg[i]], names[Fe[j]]))
            ax.plot(vel_kms[Fe[j]], col_dens[Fe[j]])
            ax.plot(vel_kms[Mg[i]], col_dens[Mg[i]])
            pdf.savefig()
            plt.close()
    #here we will try and create a scatter plot for each Mg value with the 5 Fe values (still ratio/wavelength)
    # counter = 0
    # for i in range(0, len(Mg)):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(1,1,1)
    #     plt.title(names[Mg[i]])
    #     plt.xlabel('Rest Wavelength')
    #     plt.ylabel('Ratio')
    #     for j in range(0, len(Fe)):
    #         ax.plot(lines[Fe[j]], Ratio[counter][0], 'ro')
    #         counter += 1
    #     pdf.savefig()
    #     plt.close()
os.system("open %s &" % filename)
