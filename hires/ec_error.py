from scipy.interpolate import spline
import math
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch

def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
    return cor_col

# define the data directory
dir = os.environ['HIRESDIR']

# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853
mgii2796 = 2796.3542699
#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
# array of names, lines of interest, oscillator strength:
names = ['Mg II 2796', 'Mg II 2803']
lines = [mgii2796, mgii2803]
fosc = [f2796, f2803]

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']

# define velocity ranges to plot the profile
# vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']
# define velocity ranges to plot the profile


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
            sigma[i] = flux[i] * math.sqrt(var[i])/fx[i]
            
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


## Code to combine 2796 and 2803 absorption profiles

        # 2796 on left (blueshifted side), 2803 on right (redshifted side)
        g2796 = (vel_kms[0] > -3000) & (vel_kms[0] < vflip[h])
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < 500)

        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')

        # g2796_coldens = (vel_2796 > -3000) & (vel_2796 < vflip[h])
        # g2803_coldens (vel_2803 > vflip[h]) & (vel_2803 < 500)


        

                                               ## PLOTTING BEGINS ##
    
        fig = plt.figure()

##Flux Plot
        
        ax = fig.add_subplot(2,1,1)
        
        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.legend(loc = 4)
        plt.title("Flux vs. Velocity in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        
        #Adds error bars to Flux Plot
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


##Tau Plots

        ax = fig.add_subplot(2,1,2)

        ## Need to figure out how to calculate error in Tau

        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], tau[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], tau[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        # y-axis upper lim set to 5 because no visible difference between tau = 5 and tau = infinity
        ax.set_ylim(-2, 5)
        plt.title("Tau vs. Velocity in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Tau")
        plt.legend(loc = 1)

        fig.tight_layout()
        pdf.savefig()
        plt.close()


        

## Flux Plot again, to compare to Column Density

        fig = plt.figure()
        ax = fig.add_subplot(2,1,1)

        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("Flux vs. Velocity in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Continuum Normalized Flux")
        plt.legend(loc = 4)
        
        #Adds error bars to Flux Plot
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


        
## Column Density Plots

#Try to clean up
        ax = fig.add_subplot(2,1,2)

        ax.plot(vel_2796, col_2796, color = '#2CA14B', label = names[0])
        ax.plot(vel_2803, col_2803, color = '#2C6EA1', label = names[1])
        plt.title("Column Density vs. Velocity in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Column Density (Particle/cm^2)")
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2796, col_2796, yerr = sigma_coldens2796, linewidth = 0.1, color = '#99ccff', label = 'error')
        plt.errorbar(vel_2803, col_2803, yerr = sigma_coldens2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        fig.tight_layout()
        pdf.savefig()
        plt.close()

        
os.system("open  %s &" % filename)
