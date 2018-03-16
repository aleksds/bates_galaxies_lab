#Kwamae Delva
#Makes code that plots only column density for Mg and Fe ions in order to assume N for abslin.attrib['N']
#March 9, 2018

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

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853
mgii2796 = 2796.3542699
mgi2852 = 2852.96328

feii2600 = 2600.1724835
feii2586 = 2586.6495659
feii2382 = 2382.7641781
feii2374 = 2374.4603294
feii2344 = 2344.2129601

#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
f2852 = 1.83

f2600 = 0.2394
f2586 = 0.069126
f2382 = 0.320
f2374 = 0.0313
f2344 = 0.1142

# array of names, lines of interest, oscillator strength:
lines = [mgii2796, mgii2803, mgi2852, feii2586, feii2600, feii2374, feii2382, feii2344]
names = ['MgII 2796', 'MgII 2803', 'MgI 2852', 'Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']
fosc = [f2796, f2803, f2852, f2586, f2600, f2374, f2382, f2344]

hue = ['#2CA14B', '#99ccff', '#947e94']
# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']

# define velocity ranges to plot the profile
# vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']


# define velocity ranges to plot the profile
def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
    return cor_col



                                                      ### CODE FOR AXIS AND TITLE FONT ###
                                                      
title_font = {'fontname':'Arial', 'size':'16'}
axis_font = {'fontname':'Arial', 'size':'14'}




minorLocator = AutoMinorLocator()
filename = 'Mg_Tau_Flux_Column_Comparison.pdf'
with PdfPages(filename) as pdf:
    # for h in range(0, len(gal)):
    for h in range(0, 1):
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

## Need to figure out how to calculate error in Tau
 
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


## Code to combine 2796 and 2803 absorption profiles

#FLUX Limits
        # 2796 on left (blueshifted side), 2803 on right (redshifted side)
        g2796 = (vel_kms[0] > -3000) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < 500)
        g2852 = (vel_kms[2] > -3000) & (vel_kms[2] <500)

        g2586 = (vel_kms[0] > -3000) & (vel_kms[0] <500)
        g2600 = (vel_kms[1] > -3000) & (vel_kms[1] <500)
        g2374 = (vel_kms[2] > -3000) & (vel_kms[2] <500)
        g2382 = (vel_kms[3] > -3000) & (vel_kms[3] <500)        
        g2344 = (vel_kms[4] > -3000) & (vel_kms[4] <500)

          

#COLUMN Info        
        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2852 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        col_2586 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2600 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2374 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))
        col_2382 = column(vel_kms[3],tau/(2.654E-15*fosc[3]**2 *(wave/(1+zem[h]))))
        col_2344 = column(vel_kms[4],tau/(2.654E-15*fosc[4]**2 *(wave/(1+zem[h]))))

        

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')
        vel_2852 = np.linspace(-3000,500, num = len(col_2852), endpoint = 'True')

        vel_2586 = np.linspace(-3000,500, num = len(col_2586), endpoint = 'True')
        vel_2600 = np.linspace(-3000,500, num = len(col_2600), endpoint = 'True')
        vel_2374 = np.linspace(-3000,500, num = len(col_2374), endpoint = 'True')
        vel_2382 = np.linspace(-3000,500, num = len(col_2382), endpoint = 'True')
        vel_2344 = np.linspace(-3000,500, num = len(col_2344), endpoint = 'True')
    
#Error Limits        
        sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_coldens2852 = column(vel_kms[2], sigma_coldens/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        sigma_coldens2586 = column(vel_kms[3], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2600 = column(vel_kms[4], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2374 = column(vel_kms[5], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2382 = column(vel_kms[6], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2344 = column(vel_kms[7], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        

        print('This is max coldens for MgII2796', np.array(np.max(sigma_coldens2796)))
        print('This is max coldens for MgII2803', np.array(np.max(sigma_coldens2803)))
        print('This is max coldens for MgI2852', np.array(np.max(sigma_coldens2852)))

        print('This is max coldens for MgII2586', np.array(np.max(sigma_coldens2586)))
        print('This is max coldens for MgII2600', np.array(np.max(sigma_coldens2600)))
        print('This is max coldens for MgII2374', np.array(np.max(sigma_coldens2374)))
        print('This is max coldens for MgII2382', np.array(np.max(sigma_coldens2382)))
        print('This is max coldens for MgII2344', np.array(np.max(sigma_coldens2344)))

        


os.system("open  %s &" % filename)


# To continue, modify x and y axis limits for each galaxy and wavelength and I should be good to go with calculating good N values for linetools code. Hopefully Aleks can help me out with this. 
