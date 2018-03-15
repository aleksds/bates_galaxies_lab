#Kwamae Delva  03/14/18
#Code to figure out b parameter for linetools code
#Based on kd_all_Mg_all_plots.py


import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting
from scipy.interpolate import UnivariateSpline
import pylab as pl
    
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
    for h in range(0,1):
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

        

#V_FWHM Code   from   https://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak

        # create a spline of vel_kms and flux+np.min(flux)/2 
        spline = UnivariateSpline(vel_kms[0], flux[g2796]-np.max(flux[g2796])/2, s=0)
        p1, p2 = spline.roots() # find the roots

        spline = UnivariateSpline(vel_kms[1], flux[g2803]-np.max(flux[g2803])/2, s=0)
        q1, q2 = spline.roots() # find the roots
        
        spline = UnivariateSpline(vel_kms[2], flux[g2852]-np.max(flux[g2852])/2, s=0)
        r1_, r2 = spline.roots() # find the roots

        
        pl.plot(vel_kms[0], flux[g2796])
        pl.axvspan(p1, p2, facecolor='g', alpha=0.5)
        pl.show()

        pl.plot(vel_kms[1], flux[g2803])
        pl.axvspan(q1, q2, facecolor='g', alpha=0.5)
        pl.show()

        pl.plot(vel_kms[2], flux[g2852])
        pl.axvspan(r1, r2, facecolor='g', alpha=0.5)
        pl.show()



        #Original Code from Website

        # import numpy as np
        # from scipy.interpolate import UnivariateSpline

        # def make_norm_dist(x, mean, sd):
        #     return 1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x - mean)**2/(2*sd**2))

        # x = np.linspace(10, 110, 1000)
        # green = make_norm_dist(x, 50, 10)
        # pink = make_norm_dist(x, 60, 10)

        # blue = green + pink   

        # # create a spline of x and blue-np.max(blue)/2 
        # spline = UnivariateSpline(x, blue-np.max(blue)/2, s=0)
        # r1, r2 = spline.roots() # find the roots

        # import pylab as pl
        # pl.plot(x, blue)
        # pl.axvspan(r1, r2, facecolor='g', alpha=0.5)
        # pl.show()
        

                                               ## PLOTTING BEGINS ##
    
                                               



                                        ## Combined Mg II Ion Plots##
                                        

# ##Flux Plot
        
#         fig = plt.figure()

#         ax = fig.add_subplot(2,1,1)
        
#         #Actual Flux vs. Velocity Plot
#         ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
#         ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
#         ax.set_xlim(-1500, -1000)
#         ax.set_ylim(0, 2)
#         # plt.legend(loc = 3)
#         plt.title("MgII Flux Plot for Galaxy %s" %(gal[h]), **title_font)
#         plt.xlabel("Velocity(km/s)", **axis_font)
#         plt.ylabel("C.N. Flux", **axis_font)
#         plt.rc('xtick', labelsize=10) 
#         plt.rc('ytick', labelsize=10)
        
#         #Adds error bars to Flux Plot
#         plt.errorbar(vel_kms[0][g2796], flux[g2796], yerr = sigma[g2796], label = 'error', color = '#99ccff', markevery = 10, linewidth = .1)
#         plt.errorbar(vel_kms[1][g2803], flux[g2803], yerr = sigma[g2803], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)


        
#                                                # MgI 2852 Plots ##
# # Flux Plot

#         ax = fig.add_subplot(2,1,2)
#         #Actual Flux vs. Velocity Plot
#         ax.plot(vel_kms[2], flux, linewidth=1, label = names[2], color = '#947e94')
#         ax.set_xlim(-1500, -1000)
#         ax.set_ylim(0, 2)
#         # plt.legend(loc = 3)
#         plt.title("MgI 2852 Flux Plot for Galaxy %s" %(gal[h]), **title_font)
#         plt.xlabel("Velocity(km/s)", **axis_font)
#         plt.ylabel("C.N. Flux", **axis_font)
#         plt.rc('xtick', labelsize=10) 
#         plt.rc('ytick', labelsize=10)
        
#         #Adds error bars to Flux Plot
#         # plt.errorbar(vel_kms[2][g2852], flux[g2852], yerr = sigma[g2852], label = '_nolegend_',  color = '#99ccff', markevery = 10, linewidth = .1)

#         pdf.savefig()
#         plt.close()

# os.system("open  %s &" % filename)
