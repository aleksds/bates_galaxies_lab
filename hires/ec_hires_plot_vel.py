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


# speed of light in Angstroms per second
c = const.c.to('AA/s')
c_kms = const.c.to('km/s')
#mass of electron
mass_e = 9.11**(-28)

# SDSS DATA!!!


#data directory for SDSS data
dir_SDSS = os.environ['SDSSDIR']
file = dir_SDSS + 'spec-0761-54524-0409.fits'
hdulist = fits.open(file)

coeff0_SDSS = hdulist[0].header['COEFF0']
coeff1_SDSS = hdulist[0].header['COEFF1']

flux_SDSS = hdulist[1].data

npix = len(flux_SDSS)
index = np.arange(npix)
wave_SDSS = np.array([10.**(coeff0_SDSS + coeff1_SDSS*index)])


## SDSS PLOT WORK!!


#redshift coefficient
z = 0.60329
#wavelengths:
w_int_SDSS = ['MgII2796','MgII2803', 'MgI2852', 'FeII2600']
w_rest_SDSS = np.array([2796, 2803 , 2852, 2600])
w_obs_SDSS = []
#define w_obs:
for i in range(0, len(w_rest_SDSS)):
    temp = w_rest_SDSS[i]*(1+z)
    w_obs_SDSS.append(temp)
print(w_obs_SDSS)
#wavelength to velocity conversion:
vel_obs_SDSS = np.array([])
vel_aas_SDSS = np.zeros([len(w_obs_SDSS),3842])
vel_aas_0_SDSS = []
vel_aas_1_SDSS = []
vel_aas_2_SDSS = []
vel_aas_3_SDSS = []

#vel_aas = [[vel_aas_0], [vel_aas_1], [vel_aas_2], [vel_aas_3]]
vel_kms_SDSS = np.zeros([len(w_obs_SDSS),3842])
# define the velocity scale [AA / s]
for i in range(0, len(w_obs_SDSS)):
    vel_aas_SDSS[i] = ((wave_SDSS[0]-w_obs_SDSS[i])/w_obs_SDSS[i]) * c
    # convert to [km / s]
    #vel_kms[i] = np.array([vel_aas[i].to('km/s')])
    vel_kms_SDSS[i] = vel_aas_SDSS[i]


for i in range(0, len(w_obs_SDSS)):
    vel_temp_SDSS = (w_obs_SDSS[i]) / (w_obs_SDSS[i]*(1+z)) * c
    vel_kms_temp_SDSS = np.array([vel_temp_SDSS.to('km/s')], dtype = object)
    vel_obs_SDSS = np.append(vel_obs_SDSS, vel_kms_temp_SDSS)


#Question: how long should tau be?  1 value for each line, same # as flux for each line, just the same # as flux in total?
tau_SDSS = np.zeros([len(w_obs_SDSS),len(flux_SDSS)])
# loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
for j in range(0, len(w_obs_SDSS)):
    boom = np.log(1/flux_SDSS)
    tau_SDSS[j] = boom



#  NEED TO FIND THE SDSS VERSION OF FOSC (OSCILLATOR STRENGTH) TO FINISH SDSS COL_DENS CALCULATION


    
#Column Density Calc
#col_dens_SDSS[i] = tau_SDSS[i] / (2.654E-15 * fosc[i] * (wave_SDSS/(1+z)) * fosc[i])

    
# SDSS PLOTTING

    
minorLocator = AutoMinorLocator()
filename = 'J0826_sdss.pdf'
with PdfPages(filename) as pdf:
    #prints the first graph 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(wave_SDSS[0], flux_SDSS)
    for i in range(0, len(w_rest)):
        plt.axvline(x=w_obs_SDSS[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs_SDSS[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs_SDSS[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
    plt.xlabel('Wavelength_SDSS')
    plt.ylabel('Flux_SDSS')
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
        plt.title(w_int_SDSS[i])
        ax.set_xlim(w_obs_SDSS[i]-200, w_obs_SDSS[i]+200)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.plot(wave_SDSS[0], flux_SDSS)
        plt.xlabel('Wavelength_SDSS')
        plt.ylabel('Flux_SDSS')
        plt.axvline(x=w_obs_SDSS[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=w_obs_SDSS[i]+5., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=w_obs_SDSS[i]-5., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        pdf.savefig()
        plt.close()

    #prints the graph with velocity
    for i in range(0, len(w_int)):
        fig3 = plt.figure()
        ax = fig3.add_subplot(1,1,1)
        ax.plot(vel_kms_SDSS[i], flux_SDSS)
        plt.title(w_int_SDSS[i])
        plt.axvline(x=-1600, ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=-1600+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=-1600-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.xlabel('Velocity_SDSS (Observed)')
        plt.ylabel('Flux_SDSS')
        ax.set_ylim(0., 15.)
        ax.set_xlim(-3000,500)
        pdf.savefig()
        plt.close()


# HIRES DATA!!!!


# data directory for HIRES
dir_HIRES = os.environ['HIRESDIR']

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir_HIRES+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']

for i in range(0, 1):#len(gal)):
    # read in the spectrum
    datafile = dir_HIRES+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
    data = ascii.read(datafile)
    wave_HIRES = data['wv'] * u.AA
    flux_HIRES = data['norm']  #intensity of the light you obseved   ;  normalized means you check observed light from expected light


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


##  HIRES WORK


vel_kms_HIRES = np.zeros([len(lines),len(wave_HIRES)])
# define the velocity scale [km / s]
for i in range(0, len(lines)):
    vel_kms_HIRES[i] = ((wave_HIRES-lines[i]*(1+zem[0]))/(lines[i]*(1+zem[0]))) * c_kms

#Question: how long should tau be?  1 value for each line, same # as flux for each line, just the same # as flux in total?
tau_HIRES = np.zeros([len(lines),len(flux_HIRES)])
# loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
for j in range(0, len(lines)):
    blah = np.log(1/flux_HIRES)
    tau_HIRES[j] = blah


#Column Density Calc
col_dens_HIRES[i] = tau_HIRES[i] / (2.654E-15 * fosc[i] * (wave_HIRES/(1+zem[0])) * fosc[i])

        
os.system("open %s &" % filename)
    
