#Kwamae Delva  07/27/17
#Makes Seperate Flux Plots for each Ion and compares it to column density and Voigt Profile
#Makes Combined Flux Plot and compares it to Tau and Column Density
#Code on bottom is an attempt to create functions for the plots
#Final Code

#07/31/17  =  attempt to add MgI 2852 wavelength
#08/2/17 = Mg Omnipotent Works for Mg Ions


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

        column_densities = [col_2796, col_2852, col_2852]

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')
        vel_2852 = np.linspace(-3000,500, num = len(col_2852), endpoint = 'True')

        column_velocities = [vel_2796, vel_2803, vel_2852]
        
#Error Limits        
        sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_coldens2852 = column(vel_kms[2], sigma_coldens/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))
        


                                               ## PLOTTING BEGINS ##
    
                                               

                                            ## Seperated Mg II Ion Plots##

                
        ## Mg II 2796 

        fig = plt.figure()

#Flux and Voigt Spectrum

        ax = fig.add_subplot(3,1,1)
        ax.plot(vel_kms[0], flux, linewidth=1, label = names[0], color = '#2CA14B')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("Mg II 2796 Voigt Profile for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)

        

        f = interp1d(vel_kms[0], flux)
        #Tells code to make voigt profile in velocity range -3000 to 500 when flux is below .5
        test_2796 = (vel_kms[0] > -3000) & (vel_kms[0] < 500) & (flux < 0.5)
        vel_median = np.median(vel_kms[0][test_2796])
        vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)


        boom = len(vel_kms[0][test_2796])
        one = vel_kms[0][test_2796][round(boom*0.2)]
        two = vel_kms[0][test_2796][round(boom*0.4)]
        thr = vel_kms[0][test_2796][round(boom*0.6)]
        fou = vel_kms[0][test_2796][round(boom*0.8)]

        flux_king = f(vel_new)

        xarr = vel_new
        yarr = flux_king - 1.

        voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)

                   ## Write function that combines cover fraction code with Voigt profile code ??????????? ##
                   ##Correlation between amplitude and cover frac could be key##
                   
        fitter = fitting.LevMarLSQFitter()
        voi_fit = fitter(voi_init, xarr, yarr)

        ax.plot(xarr,voi_fit(xarr)+1, color='red')
        plt.title("Mg 2796 Voigt Profile for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim (0,2)
        ax.set_xlim(-3000,500)

## Column Density Plot

        ax = fig.add_subplot(3,1,3)
        ax.plot(vel_2796, col_2796, linewidth =1, color = '#2CA14B', label = names[0])
        plt.title("MgII 2796 Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)",**axis_font)
        plt.ylabel("Column Density", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        # plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2796, col_2796, yerr = sigma_coldens2796, linewidth = 0.1, color = '#99ccff', label = 'error')
        
        # fig.tight_layout()
        pdf.savefig()
        plt.close()



        
#### Mg II 2803

        fig = plt.figure()

# Flux and Voigt Profile
        
        ax = fig.add_subplot(3,1,1)
        ax.plot(vel_kms[1], flux, linewidth=1, label = names[0], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("Mg II 2803 Voigt Profile for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)

        f = interp1d(vel_kms[1], flux)
        test_2803 = (vel_kms[1] > -3000) & (vel_kms[1] < 500) & (flux < 0.5)
        vel_median = np.median(vel_kms[1][test_2803])
        vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)


        boom = len(vel_kms[1][test_2803])
        one = vel_kms[1][test_2803][round(boom*0.2)]
        two = vel_kms[1][test_2803][round(boom*0.4)]
        thr = vel_kms[1][test_2803][round(boom*0.6)]
        fou = vel_kms[1][test_2803][round(boom*0.8)]

        flux_king = f(vel_new)

        xarr = vel_new
        yarr = flux_king - 1.

        voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)
        fitter = fitting.LevMarLSQFitter()
        voi_fit = fitter(voi_init, xarr, yarr)

        ax.plot(xarr,voi_fit(xarr)+1, color='red')
        plt.title("Mg 2803 Voigt Profile for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux",**axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim (0,2)
        ax.set_xlim(-3000,500)
        
## Column Density Plot
        
        ax = fig.add_subplot(3,1,3)
        ax.plot(vel_2803, col_2803, linewidth =1, color = '#2C6EA1', label = names[1])
        plt.title("MgII 2803 Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)",**axis_font)
        plt.ylabel("Column Density",**axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        # plt.legend(loc = 1)

        #Adds error bars to plots
        plt.errorbar(vel_2803, col_2803, yerr = sigma_coldens2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        # fig.tight_layout()
        pdf.savefig()
        plt.close()


                                        ## Combined Mg II Ion Plots##
                                        

##Flux Plot
        
        fig = plt.figure()

        ax = fig.add_subplot(3,1,1)
        
        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        # plt.legend(loc = 3)
        plt.title("MgII Flux Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        
        #Adds error bars to Flux Plot
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


##Tau Plot        
        ax = fig.add_subplot(3,1,3)

        ## Need to figure out how to calculate error in Tau

        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], tau[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], tau[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        # y-axis upper lim set to 5 because no visible difference between tau = 5 and tau = infinity
        ax.set_ylim(-.2, 5)
        plt.title("MgII Tau Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("Tau", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        # plt.legend(loc = 4)

        # fig.tight_layout()
        pdf.savefig()
        plt.close()


        

## Flux Plot again, to compare to Column Density

##Could possible change this to Voigt profile plot of combined MgII ions flux plot
    
        fig = plt.figure()

        ax = fig.add_subplot(3,1,1)

        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("MgII Flux Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        # plt.legend(loc = 3)
        
        #Adds error bars to Flux Plot
        plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
        plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


        
## Column Density Plots

#Try to clean up by removing extra peaks beside 1250 km/s

        ax = fig.add_subplot(3,1,3)

        ax.plot(vel_2796, col_2796, color = '#2CA14B', label = names[0])
        ax.plot(vel_2803, col_2803, color = '#2C6EA1', label = names[1])
        plt.title("MgII Col. Dens. Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity (km/s)", **axis_font)
        plt.ylabel("Column Density", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        ax.set_ylim(0, 5E12)
        ax.set_xlim(-3000,500)
        # plt.legend(loc = 1)
        
        #Adds error bars to plots
        plt.errorbar(vel_2796, col_2796, yerr = sigma_coldens2796, linewidth = 0.1, color = '#99ccff', label = 'error')
        plt.errorbar(vel_2803, col_2803, yerr = sigma_coldens2803, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        # fig.tight_layout()
        pdf.savefig()
        plt.close()

        
                                               # MgI 2852 Plots ##
# Flux Plot

        fig = plt.figure()

        ax = fig.add_subplot(3,1,1)
        #Actual Flux vs. Velocity Plot
        ax.plot(vel_kms[2], flux, linewidth=1, label = names[2], color = '#947e94')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        # plt.legend(loc = 3)
        plt.title("MgI 2852 Flux Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("C.N. Flux", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        
        #Adds error bars to Flux Plot
        # plt.errorbar(vel_kms[2][g2852], flux[g2852], yerr = sigma[g2852], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


# Tau Plot

        ax = fig.add_subplot(3,1,3)
        ax.plot(vel_kms[2], tau, linewidth=1, label = names[2], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        # y-axis upper lim set to 5 because no visible difference between tau = 5 and tau = infinity
        ax.set_ylim(-.2, 5)
        plt.title("MgI 2852 Tau Plot for Galaxy %s" %(gal[h]), **title_font)
        plt.xlabel("Velocity(km/s)", **axis_font)
        plt.ylabel("Tau", **axis_font)
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
        # plt.legend(loc = 4)


# No Voigt Profile Because Flux plot doesn't have significant absorption line for MgI 2852
        
## No Column Density Plots because there's barely any MgI 2852 Present


        # ax = fig.add_subplot(3,1,3)
        # ax.plot(vel_2852, col_2852, color = '#2C6EA1', label = names[1])
        # plt.title("MgI 2852 Column Density Plot for Galaxy %s" %(gal[h]), **title_font)
        # plt.xlabel("Velocity (km/s)", **axis_font)
        # plt.ylabel("Column Density (Particle/cm^2)", **axis_font)
        # plt.rc('xtick', labelsize=10) 
        # plt.rc('ytick', labelsize=10)
        # ax.set_ylim(0, 5E12)
        # ax.set_xlim(-3000,500)
        # # plt.legend(loc = 1)
        
        # #Adds error bars to plots
        # # plt.errorbar(vel_2852, col_2852, yerr = sigma_coldens2852, linewidth = 0.1, color = '#99ccff', label = '_nolegend_')

        # fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("open  %s &" % filename)















                                                              # ##### Work In Progress####



# ##Add limit of range to below functions (g2796 etc)

# def Voigt(velocityy, fluxx, namess):

#     fig = plt.figure
    
#     for i in len(namess):
        
#         ax = fig.add_subplot(2,1,2)
        
#         f = interp1d(velocityy[i], fluxx)
#         if fluxx < .5:
            
#             #Tells code to make Voigt Profile in velocity range -3000 to 500 when flux is below .5
#             to_fit = (velocityy[i] > -3000) & (velocityy[i] < 500) & (fluxx < .5)

#             #finds median of velocities that meet to_fit parameter
#             vel_median = np.median(velocityy[i][to_fit])

#             #insertes evenly spaced integers in range 1000 less/above vel_median
#             vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)

#             #seperates velocities into 4 even quadrants ; round just rounds the number to an integer
#             boom = len(velocityy[i][to_fit])
#             one = velocityy[i][to_fit][round(boom*0.2)]
#             two = velocityy[i][to_fit][round(boom*0.4)]
#             thr = velocityy[i][to_fit][round(boom*0.6)]
#             fou = velocityy[i][to_fit][round(boom*0.8)]

#             #Voigt Profile Parameter Info
#             flux_king = f(vel_new)
#             xarr = vel_new
#             yarr = flux_king -1

#             #Makes Voigt Profile with amplitude, center, and fwhm_G/L as parameters
#             voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)
#             fitter = fitting.LevMarLSQFitter()
#             voi_fit = fitter(voi_init, xarr, yarr)

#             #Plots Voigt Profile in red
#             ax.plot(xarr,voi_fit(xarr)+1, color='red')
#             plt.xlabel("Velocity(km/s)")
#             plt.ylabel("Continuum Normalized Flux")
#             ax.set_ylim (0,2)
#             ax.set_xlim(-3000,500)

#             fig.tight_layout()
#             pdf.savefig()
#             plt.close()

#         else:
#             print( "No Voigt Profile Possible for %s," %(names[i]), "in Galaxy %s" %(gal[h]))

# def COLUMN(velocityy, col_denss, namess, hue):
#     fig = plt.figure()
#     for i in len(namess):
#         ax = fig.add_subplot(2,1,2)
#         ax.plot(velocityy[i], col_dens[i], linewidth =1, color = hue[i], label = namess[i])
#         plt.title("Column Density Plot for Galaxy %s" %(gal[h]))
#         plt.xlabel("Velocity (km/s)")
#         plt.ylabel("Column Density (Particle/cm^2)")
#         ax.set_ylim(0, 5E12)
#         ax.set_xlim(-3000,500)
#         plt.legend(loc = 1)

#         fig.tight_layout()
#         pdf.savefig()
#         plt.close()

# def FLUX(velocityy, fluxx, namess, hue):
#     fig = plt.figure()
#     for i in len(namess):

#         ax = fig.add_subplot(2,1,1)
#         ax.plot(velocityy[i], fluxx, linewidth=1, label = namess[i], color = hue[i])
#         ax.set_xlim(-3000, 500)
#         ax.set_ylim(0, 2)
#         plt.legend(loc = 3)
#         plt.title("Flux Plot for Galaxy %s" %(gal[h]))
#         plt.xlabel("Velocity(km/s)")
#         plt.ylabel("Continuum Normalized Flux")

#         fig.tight_layout()
#         pdf.savefig()
#         plt.close()

# def TAU(velocityy, tauu, namess, hue, limit):

#     fig = plt.figure()

#     for i in len(namess):

#         ax = fig.add_subplot(2,1,2)

#         #Actual Flux vs. Velocity Plot
#         ax.plot(velocityy[i][limit[i]], tau[i][limit[i]], linewidth=1, label = namess[i], color = hue[i])
#         ax.set_xlim(-3000, 500)
#         # y-axis upper lim set to 5 because no visible difference between tau = 5 and tau = infinity
#         ax.set_ylim(-.2, 5)
#         plt.title("Tau Plot for Galaxy %s" %(gal[h]))
#         plt.xlabel("Velocity(km/s)")
#         plt.ylabel("Tau")
#         plt.legend(loc = 4)
        
#         fig.tight_layout()
#         pdf.savefig()
#         plt.close()

# def COMB_FLUX(velocityy, fluxx, namess, hue, limit):
#     fig = plt.figure()

#     ax = fig.add_subplot(2,1,1)

#     #Actual Flux vs. Velocity Plot
#     ax.plot(velocityy[0][limit[0]], flux[limit[0]], linewidth=1, label = namess[0], color = hue[0])
#     ax.plot(velocityy[1][limit[1]], flux[limit[1]], linewidth=1, label = namess[1], color = hue[1])
#     ax.set_xlim(-3000, 500)
#     ax.set_ylim(0, 2)
#     plt.title("MgII Doublet Flux Plot for Galaxy %s" %(gal[h]))
#     plt.xlabel("Velocity(km/s)")
#     plt.ylabel("Continuum Normalized Flux")
#     plt.legend(loc = 3)
    
#     fig.tight_layout()
#     pdf.savefig()
#     plt.close()

# #Calling Functions
#     # FLUX(vel_kms, flux, names, hue)
#     # TAU(vel_kms, tau, names, hue, limit)
#     # COLUMN(column_velocities, column_densities, names, hue)
#     # Voigt(vel_kms, flux, names)
#     # COMB_FLUX(vel_kms, flux, names, hue, limit)
