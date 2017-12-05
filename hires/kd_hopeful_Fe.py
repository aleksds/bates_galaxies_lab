#Kwamo 07/31/17
#All plots for all the Fe ions, grouped together for comparison
#08/03/17 Iron Code Works

#Future Steps: add parameter for optical thickness 

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
# hue = ['b8860b', '708090', '228b22', 'ffe4c4', '8a2be2']
hue = ['green', 'magenta', 'blue', 'brown', 'black']

# arrays of galaxy names, redshifts, approximate centroid velocities
gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']

# define velocity ranges to plot the profile
# vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']

def column(vel, col_dens):
    cor_col = np.array([])
    for i in range(0, len(vel)):
        if (vel[i] >= -3000 and vel[i] <= 500):
            cor_col = np.append(cor_col, col_dens[i])
    return cor_col





title_font = {'fontname':'Arial', 'size':'10'}
axis_font = {'fontname':'Arial', 'size':'8'}




## Function for Plotting
# def Omnipotent(

minorLocator = AutoMinorLocator()
filename = 'Fe_Tau_Flux_Coldens_Comparison.pdf'
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

#Velocity Info

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5

## Tau Info

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        blah = np.log(1/flux)
        tau = blah

        # tau_limiter = np.array([])
        # for i in vel_kms:
        #     limit_calc = (vel_kms[i] > -3000) & (vel_kms[i] <500)
        #     tau_limiter = np.append(tau_limiter, limit_calc)

            
        g2586 = (vel_kms[0] > -3000) & (vel_kms[0] <500)
        g2600 = (vel_kms[1] > -3000) & (vel_kms[1] <500)
        g2374 = (vel_kms[2] > -3000) & (vel_kms[2] <500)
        g2382 = (vel_kms[3] > -3000) & (vel_kms[3] <500)        
        g2344 = (vel_kms[4] > -3000) & (vel_kms[4] <500)

        tau_limit = [g2586, g2600, g2374, g2382, g2344]
        
#Column Density Info

        ## Column Densities
        # column_dens = np.array([])
        # for i in range(0, len(names)):
        #     column_calc = column(vel_kms[i],tau/(2.654E-15*fosc[i]**2 *(wave/(1+zem[h]))))
        #     column_dens = np.append(column_dens, column_calc)

            
        col_2586 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2600 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2374 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))
        col_2382 = column(vel_kms[3],tau/(2.654E-15*fosc[3]**2 *(wave/(1+zem[h]))))
        col_2344 = column(vel_kms[4],tau/(2.654E-15*fosc[4]**2 *(wave/(1+zem[h]))))
        
        column_densities = [col_2586, col_2600,col_2374, col_2382, col_2344]

        ## Velocity for Column Density
        # column_velo = np.array ([])
        # for i in column_dens:
        #     colvel_calc = np.linspace(-3000,500, num = len([i]), endpoint = 'True')
        #     column_velo = np.append(column_velo, colvel_calc)
            
        vel_2586 = np.linspace(-3000,500, num = len(col_2586), endpoint = 'True')
        vel_2600 = np.linspace(-3000,500, num = len(col_2600), endpoint = 'True')
        vel_2374 = np.linspace(-3000,500, num = len(col_2374), endpoint = 'True')
        vel_2382 = np.linspace(-3000,500, num = len(col_2382), endpoint = 'True')
        vel_2344 = np.linspace(-3000,500, num = len(col_2344), endpoint = 'True')

        column_velocities = [vel_2586, vel_2600, vel_2374, vel_2382, vel_2344]

                                                 ##PLOTTING BEGINS##

#Flux and Voigt Spectrum

# # Voigt Profile

#         fig = plt.figure()
#         for i in range(0, len(names)):
#             if np.any(flux < .5):

#                 ax = fig.add_subplot(3,2,i+1)
#                 ax.plot(vel_kms[i], flux, linewidth=1, label = names[i], color = hue[i])
#                 ax.set_xlim(-3000, 500)
#                 ax.set_ylim(0, 2)
#                 plt.title(" %s Flux Plot %s" %((names[i]), (gal[h]),**title_font))
#                 plt.xlabel("Velocity(km/s)",**axis_font)
#                 plt.ylabel("C.N. Flux",**axis_font)

            
#                 f = interp1d(vel_kms[i], flux)
#                 #Tells code to make Voigt Profile in velocity range -3000 to 500 when flux is below .5
#                 to_fit = (vel_kms[i] > -3000) & (vel_kms[i] < 500) & (flux < .5)

#                 #finds median of velocities that meet to_fit parameter
#                 vel_median = np.median(vel_kms[i][to_fit])
                    
#                 #insertes evenly spaced integers in range 1000 less/above vel_median
#                 vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)
                    
#                 #seperates velocities into 4 even quadrants ; round just rounds the number to an integer
#                 boom = len(vel_kms[i][to_fit])
#                 one = vel_kms[i][to_fit][round(boom*0.2)]
#                 two = vel_kms[i][to_fit][round(boom*0.4)]
#                 thr = vel_kms[i][to_fit][round(boom*0.6)]
#                 fou = vel_kms[i][to_fit][round(boom*0.8)]

#                 #Voigt Profile Parameter Info
#                 flux_king = f(vel_new)
#                 xarr = vel_new
#                 yarr = flux_king -1
                
#                 #Makes Voigt Profile with amplitude, center, and fwhm_G/L as parameters
#                 voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)
#                 fitter = fitting.LevMarLSQFitter()
#                 voi_fit = fitter(voi_init, xarr, yarr)

#                 #Plots Voigt Profile in red
#                 ax.plot(xarr,voi_fit(xarr)+1, color='red')
#                 plt.xlabel("Velocity(km/s)",**axis_font)
#                 plt.ylabel("C.N. Flux",**axis_font)
#                 ax.set_ylim (0,2)
#                 ax.set_xlim(-3000,500)

#         fig.tight_layout()
#         pdf.savefig()
#         plt.close()

#             else:
#                 fig = plt.figure()
#                 for i in range(0, len(names)):

#                     ax = fig.add_subplot(3,2,i+1)
            
#                     ax.plot(vel_kms[i], flux, linewidth=1, label = names[i], color = hue[i])
#                     ax.set_xlim(-3000, 500)
#                     ax.set_ylim(0, 2)
#                     plt.title(" %s Flux Plot %s" %((names[i]), (gal[h]),**title_font))
#                     plt.xlabel("Velocity(km/s)",**axis_font)
#                     plt.ylabel("Flux",**axis_font)
#         fig.tight_layout()
#         pdf.savefig()
#         plt.close()


#Flux, gonna remove once voigt works
        fig = plt.figure()
        plt.rc('xtick', labelsize=6) 
        plt.rc('ytick', labelsize=6)
        for i in range(0, len(names)):

            ax = fig.add_subplot(3,2,i+1)
            
            ax.plot(vel_kms[i], flux, linewidth=1, label = names[i], color = hue[i])
            ax.set_xlim(-3000, 500)
            ax.set_ylim(0, 2)
            plt.title(" %s Flux Plot %s" %((names[i]), (gal[h])), **title_font)
            plt.xlabel("Velocity(km/s)",**axis_font)
            plt.ylabel("C.N. Flux",**axis_font)
            # plt.rc('xtick', labelsize=6) 
            # plt.rc('ytick', labelsize=6)
            
        fig.tight_layout()
        pdf.savefig()
        plt.close()

        
## Column Density Plot

        fig = plt.figure()
        for i in range(0, len(names)):
            ax = fig.add_subplot(3,2,i+1)
            ax.plot(column_velocities[i], column_densities[i], linewidth =1, color = hue[i], label = names[i])
            plt.title("%s Col Dens. Plot %s" %((names[i]), (gal[h])), **title_font)
            plt.xlabel("Velocity (km/s)", **axis_font)
            plt.ylabel("Col. Dens.",**axis_font)
            ax.set_ylim(0, 2E13)
            ax.set_xlim(-3000,500)
            plt.rc('xtick', labelsize=6) 
            plt.rc('ytick', labelsize=6)
            # plt.legend(loc = 1)
            
        fig.tight_layout()
        pdf.savefig()
        plt.close()

## Tau Plot

        fig = plt.figure()
        for i in range(0, len(names)):
            
            ax = fig.add_subplot(3,2,i+1)
            ax.plot(vel_kms[i][tau_limit[i]], tau[tau_limit[i]], linewidth=1, label = names[i], color = hue[i])
            ax.set_xlim(-3000, 500)
            # y-axis upper lim set to 5 because no visible difference between tau = 5 and tau = infinity
            ax.set_ylim(-.2, 5)
            plt.title("%s Tau Plot %s" %((names[i]), (gal[h])),**title_font)
            plt.xlabel("Velocity(km/s)",**axis_font)
            plt.ylabel("Tau",**axis_font)
            plt.rc('xtick', labelsize=6) 
            plt.rc('ytick', labelsize=6)
            # plt.legend(loc = 4)

        fig.tight_layout()
        pdf.savefig()
        plt.close()
        
        
os.system("open  %s &" % filename)
