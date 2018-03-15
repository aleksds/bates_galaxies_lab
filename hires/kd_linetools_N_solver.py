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

#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
f2852 = 1.83

# array of names, lines of interest, oscillator strength:
lines = [mgii2796, mgii2803, mgi2852]
names = ['MgII 2796', 'MgII 2803', 'MgI 2852']
fosc = [f2796, f2803, f2852]

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
        
        limit = [g2796, g2803, g2852]        

#COLUMN Info        
        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2852 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        # print('This is max coldens for MgII2796', np.array(np.max(col_2796)))
        # print('This is max coldens for MgII2803', np.array(np.max(col_2803)))
        # print('This is max column density for MgI2852', np.array(np.max(col_2852))) 

        column_densities = [col_2796, col_2852, col_2852]

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')
        vel_2852 = np.linspace(-3000,500, num = len(col_2852), endpoint = 'True')

        column_velocities = [vel_2796, vel_2803, vel_2852]
        
#Error Limits        
        sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_coldens2852 = column(vel_kms[2], sigma_coldens/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        print('This is max coldens for MgII2796', np.array(np.max(sigma_coldens2796)))
        print('This is max coldens for MgII2803', np.array(np.max(sigma_coldens2803)))
        print('This is max column density for MgI2852', np.array(np.max(sigma_coldens2852))) 
        


                                               ## PLOTTING BEGINS ##

## Seperated Mg II Ion Plots##

                
        ## Mg II 2796 

        fig = plt.figure()
## Column Density Plot

        ax = fig.add_subplot(2,1,1)
        ax.plot(vel_2796, col_2796, linewidth =1, color = '#2CA14B', label = names[0])
        plt.title("MgII 2796 Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)",**axis_font)
        plt.ylabel("Column Density", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 2E14)
        ax.set_xlim(-1500,-1000)
        # plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2796, col_2796, yerr = sigma_coldens2796, linewidth = 0.1, color = '#99ccff', label = 'error')


        ##Mg II 2803
        
## Column Density Plot
        
        ax = fig.add_subplot(2,1,2)
        ax.plot(vel_2803, col_2803, linewidth =1, color = '#2C6EA1', label = names[1])
        plt.title("MgII 2803 Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)",**axis_font)
        plt.ylabel("Column Density",**axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 7.5E14)
        ax.set_xlim(-2250,-1950)
        # plt.legend(loc = 1)

        #Adds error bars to plots
        plt.errorbar(vel_2803, col_2803, yerr = sigma_coldens2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        # fig.tight_layout()
        pdf.savefig()
        plt.close()



        ## Combined Mg II Ion Plots##

## Column Density Plots

        fig=plt.figure()

        ax = fig.add_subplot(2,1,1)

        ax.plot(vel_2796, col_2796, color = '#2CA14B', label = names[0])
        ax.plot(vel_2803, col_2803, color = '#2C6EA1', label = names[1])
        plt.title("MgII Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)", **axis_font)
        plt.ylabel("Column Density", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 7.5E14)
        ax.set_xlim(-3000,500)
        # plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2796, col_2796, yerr = sigma_coldens2796, linewidth = 0.1, color = '#99ccff', label = 'error')
        plt.errorbar(vel_2803, col_2803, yerr = sigma_coldens2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')


        ##Mg I 2852

## No Column Density Plots because there's barely any MgI 2852 Present


        ax = fig.add_subplot(2,1,2)
        ax.plot(vel_2852, col_2852, color = '#2C6EA1', label = names[1])
        plt.title("MgI 2852 Column Density Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)", **axis_font)
        plt.ylabel("Column Density (Particle/cm^2)", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, .6E11)
        ax.set_xlim(-3000,500)
        # plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2852, col_2852, yerr = sigma_coldens2852, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("open  %s &" % filename)


# To continue, modify x and y axis limits for each galaxy and wavelength and I should be good to go with calculating good N values for linetools code. Hopefully Aleks can help me out with this. 
