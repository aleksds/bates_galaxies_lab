# Makes plot for covering fraction vs velocity


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
vcen = gal_info['vcen']
vmax = gal_info['vmax']
vflip = gal_info['vflip']



title_font = {'fontname':'Arial', 'size':'16'}
axis_font = {'fontname':'Arial', 'size':'14'}




minorLocator = AutoMinorLocator()
filename = 'Cover_Fraction.pdf'
with PdfPages(filename) as pdf:
    for h in range(0, len(gal)):
    # for h in range(0, 1):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5
    
        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        # plots the profiles where the 2796 profile on the blue side and the 2803 profile on the red side
        
        # g2796 = (vel_kms[0] > -3000) & (vel_kms[0] < vflip[h])
        g2796 = (vel_kms[0] > -3000) & (vel_kms[0] <500)
        g2803 = (vel_kms[1] > vflip[h]) & (vel_kms[1] < 500)
        # g2803 = (vel_kms[1] > -3000) & (vel_kms[1] < 500)

        fig = plt.figure()

# ##Flux Plots
        
#         ax = fig.add_subplot(3,1,1)
#         plt.title("Flux Plot in Galaxy %s" %(gal[h]))
#         plt.xlabel("Velocity(km/s)")
#         plt.ylabel("C.N. Flux")
#         ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
#         ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
#         ax.set_xlim(-3000, 500)
#         ax.set_ylim(0, 2)

##Voigt Profile


        ax = fig.add_subplot(3,1,1)
        ax.plot(vel_kms[0][g2796], flux[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], flux[g2803], linewidth=1, label = names[0], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)
        plt.title("Mg Dub Voigt in %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("C.N. Flux")








        ##2796 Segment

        f = interp1d(vel_kms[0], flux)
        #Tells code to make voigt profile in velocity range -3000 to 500 when flux is below .5
        test_2796 = (vel_kms[0] > -3000) & (vel_kms[0] < 500) & (flux < 0.5)
        vel_median = np.median(vel_kms[0][test_2796])
        vel_new = np.linspace(vel_median-1000, vel_median+1000, num=2001, endpoint=True)


        # boom = len(vel_kms[0][test_2796])
        # one = vel_kms[0][test_2796][round(boom*0.1)]
        # two = vel_kms[0][test_2796][round(boom*0.2)]
        # thr = vel_kms[0][test_2796][round(boom*0.3)]
        # fou = vel_kms[0][test_2796][round(boom*0.4)]
        # fiv = vel_kms[0][test_2796][round(boom*0.5)]
        # six = vel_kms[0][test_2796][round(boom*0.6)]
        # sev = vel_kms[0][test_2796][round(boom*0.7)]
        # eig = vel_kms[0][test_2796][round(boom*0.8)]
        # nin = vel_kms[0][test_2796][round(boom*0.9)]

        # flux_king = f(vel_new)

        # xarr = vel_new
        # yarr = flux_king - 1.
        

        # voi_init = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fiv, fwhm_L=fiv-fou, fwhm_G=fiv-fou)+Voigt1D(amplitude_L=-1.0, x_0=six, fwhm_L=six-fiv, fwhm_G=six-fiv)+Voigt1D(amplitude_L=-1.0, x_0=sev, fwhm_L=sev-six, fwhm_G=sev-six)+Voigt1D(amplitude_L=-1.0, x_0=eig, fwhm_L=eig-sev, fwhm_G=eig-sev)+Voigt1D(amplitude_L=-1.0, x_0=nin, fwhm_L=nin-eig, fwhm_G=nin-eig)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)
        
        #            ## Write function that combines cover fraction code with Voigt profile code ??????????? ##
        #            ##Correlation between amplitude and cover frac could be key##
                   
        # fitter = fitting.LevMarLSQFitter()
        # voi_fit = fitter(voi_init, xarr, yarr)


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

        









        ##2803 Segment

        
        f_2803 = interp1d(vel_kms[1], flux)
        test_2803 = (vel_kms[1] > -3000) & (vel_kms[1] < 500) & (flux < 0.5)
        vel_median_2803 = np.median(vel_kms[1][test_2803])
        vel_new_2803 = np.linspace(vel_median_2803-1000, vel_median_2803+1000, num=2001, endpoint=True)


## increasing number of parameters actually doesn't make graph more defined, but worse.

        # boom_2803 = len(vel_kms[1][test_2803])
        # one_2803 = vel_kms[1][test_2803][round(boom_2803*0.1)]
        # two_2803 = vel_kms[1][test_2803][round(boom_2803*0.2)]
        # thr_2803 = vel_kms[1][test_2803][round(boom_2803*0.3)]
        # fou_2803 = vel_kms[1][test_2803][round(boom_2803*0.4)]
        # fiv_2803 = vel_kms[1][test_2803][round(boom_2803*0.5)]
        # six_2803 = vel_kms[1][test_2803][round(boom_2803*0.6)]
        # sev_2803 = vel_kms[1][test_2803][round(boom_2803*0.7)]
        # eig_2803 = vel_kms[1][test_2803][round(boom_2803*0.8)]
        # nin_2803 = vel_kms[1][test_2803][round(boom_2803*0.9)]

        # flux_king_2803 = f_2803(vel_new_2803)

        # xarr_2803 = vel_new_2803
        # yarr_2803 = flux_king_2803 - 1.

        
        # voi_init_2803 = Voigt1D(amplitude_L=-1.0, x_0=one, fwhm_L=two-one, fwhm_G=two-one)+Voigt1D(amplitude_L=-1.0, x_0=two, fwhm_L=thr-two, fwhm_G=thr-two)+Voigt1D(amplitude_L=-1.0, x_0=thr, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fou, fwhm_L=fou-thr, fwhm_G=fou-thr)+Voigt1D(amplitude_L=-1.0, x_0=fiv_2803, fwhm_L=fiv_2803-fou, fwhm_G=fiv_2803-fou)+Voigt1D(amplitude_L=-1.0, x_0=six_2803, fwhm_L=six_2803-fiv_2803, fwhm_G=six_2803-fiv_2803)+Voigt1D(amplitude_L=-1.0, x_0=sev_2803, fwhm_L=sev_2803-six_2803, fwhm_G=sev_2803-six_2803)+Voigt1D(amplitude_L=-1.0, x_0=eig_2803, fwhm_L=eig_2803-sev_2803, fwhm_G=eig_2803-sev_2803)+Voigt1D(amplitude_L=-1.0, x_0=nin_2803, fwhm_L=nin_2803-eig_2803, fwhm_G=nin_2803-eig_2803)+Voigt1D(amplitude_L=-1.0, x_0=vel_median_2803, fwhm_L=200, fwhm_G=200)
        
        # fitter_2803 = fitting.LevMarLSQFitter()
        # voi_fit_2803 = fitter_2803(voi_init_2803, xarr_2803, yarr_2803)


        boom_2803 = len(vel_kms[1][test_2803])
        one_2803 = vel_kms[1][test_2803][round(boom_2803*0.2)]
        two_2803 = vel_kms[1][test_2803][round(boom_2803*0.4)]
        thr_2803 = vel_kms[1][test_2803][round(boom_2803*0.6)]
        fou_2803 = vel_kms[1][test_2803][round(boom_2803*0.8)]

        flux_king_2803 = f_2803(vel_new_2803)

        xarr_2803 = vel_new_2803
        yarr_2803 = flux_king_2803 - 1.

        voi_init_2803 = Voigt1D(amplitude_L=-1.0, x_0=one_2803, fwhm_L=two_2803-one_2803, fwhm_G=two_2803-one_2803)+Voigt1D(amplitude_L=-1.0, x_0=two_2803, fwhm_L=thr_2803-two_2803, fwhm_G=thr_2803-two_2803)+Voigt1D(amplitude_L=-1.0, x_0=thr_2803, fwhm_L=fou_2803-thr_2803, fwhm_G=fou_2803-thr_2803)+Voigt1D(amplitude_L=-1.0, x_0=fou_2803, fwhm_L=fou_2803-thr_2803, fwhm_G=fou_2803-thr_2803)+Voigt1D(amplitude_L=-1.0, x_0=vel_median, fwhm_L=200, fwhm_G=200)

        fitter_2803 = fitting.LevMarLSQFitter()
        voi_fit_2803 = fitter_2803(voi_init_2803, xarr_2803, yarr_2803)


                                                       ### LINES ON PLOT CODE ###
        
        ax.plot(xarr,voi_fit(xarr)+1, color='red')    #Red line is the voigt profile fit for MgII 2796
        ax.plot(xarr_2803,voi_fit_2803(xarr_2803)+1, color='magenta')     #Magenta line is for MgII 2803??
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("C.N. Flux")
        ax.set_ylim (0,2)
        ax.set_xlim(-3000,500)



        

## Covering Fraction Plot

        covrfrac = 1-flux

        ax = fig.add_subplot(3,1,2)
        plt.title("Covering Frac. in %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Covering Fraction")
        ax.plot(vel_kms[0][g2796], covrfrac[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], covrfrac[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        ax.set_ylim(0, 2)

## Tau Plots

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line
        #Used 2803 line because all graphs have their max absorption in Mg2803 wavelength region
        #f_new_2796 = ((flux - np.min([g2796]))/(1 - np.min(flux[g2796])))
        #f_new_2803 = ((flux - np.min(flux[g2803]))/(1 - np.min(flux[g2803])))

        f_new_2796 = ((flux - np.min(voi_fit(xarr)+1))/(1 - np.min(flux[g2796])))
        f_new_2803 = ((flux - np.min(voi_fit_2803(xarr_2803)+1))/(1 - np.min(flux[g2803])))
        
        # blah = np.log(1/flux)
        tau_2796 = np.log (1/f_new_2796)
        tau_2803 = np.log (1/f_new_2803)

        
        ax = fig.add_subplot(3,1,3)
        plt.title("Tau in Galaxy %s" %(gal[h]))
        plt.xlabel("Velocity(km/s)")
        plt.ylabel("Tau")
        ax.plot(vel_kms[0][g2796], tau_2796[g2796], linewidth=1, label = names[0], color = '#2CA14B')
        ax.plot(vel_kms[1][g2803], tau_2803[g2803], linewidth=1, label = names[1], color = '#2C6EA1')
        ax.set_xlim(-3000, 500)
        # y-axis limit set to 5 because tau of 5 and t of infinity have no visible difference
        ax.set_ylim(0, 6)


        fig.tight_layout()
        pdf.savefig()
        plt.close()
os.system("open %s &" % 'Cover_Fraction.pdf')
