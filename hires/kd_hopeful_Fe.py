#Kwamo 07/31/17
#Beginning of Fe code, based on kd_omnipotent_Mg

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting

# define the data directory
dir = os.environ['HIRESDIR']

#wavelengths of relevant absorption lines

feii2600 = 2600.1724835
feii2586 = 2586.6495659
feii2382 = 2382.7641781
feii2374 = 2374.4603294
feii2344 = 2344.2129601

#oscillator strengths

f2600 = 0.2394
f2586 = 0.069126
f2382 = 0.320
f2374 = 0.0313
f2344 = 0.1142

#Code Simplifiers
lines = [feii2586, feii2600, feii2374, feii2382, feii2344]
names = ['Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']
fosc = [f2586, f2600, f2374, f2382, f2344]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']

# define velocity ranges to plot the profile
# vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']
# define velocity ranges to plot the profile

## Function for Plotting
def Omnipotent(

minorLocator = AutoMinorLocator()
filename = 'Tau_Flux_Coldens_Comparison.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']
        
## Error Bar Calculations

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * np.sqrt(var[i])/fx[i]
            
#Error in Column Density Calculations

        sigma_coldens = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma_coldens[i] = (sigma[i]/flux[i])

#Velocity Calculations

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5

## Tau Calculations

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        blah = np.log(1/flux)
        tau = blah

