# Kwamae Delva
# Testing linetools example functions

import astropy.units as u
from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio
from linetools.analysis import voigt as lav
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectra.xspectrum1d import XSpectrum1D
import imp
import warnings
from astropy.io import ascii
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting

# plots
from matplotlib import pyplot as plt
import pylab
# pylab.rcParams['figure.figsize'] = (8.0, 6.0)
import matplotlib
matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10)


# Define a plotting function
def plt_line(spec):
    plt.clf()
    plt.figure(dpi=1200)
    plt.plot(spec.dispersion.value, spec.flux.value, 'k-', drawstyle='steps-mid', lw=1.5)
    plt.xlim(3644., 3650.)
    plt.ylabel('Normalized Flux', fontsize=20.)
    plt.xlabel('Wavelength', fontsize=20.)
    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.MultipleLocator(2.))
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.ylim(0., 1.1)
    plt.show()
    plt.close()


try:
    import seaborn as sns; sns.set(context="notebook",font_scale=2)
except:
    pass

warnings.filterwarnings('ignore')

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
filename = 'Linetools_Voigt_and_Col_Dens.pdf'
with PdfPages(filename) as pdf:
    # for h in range(0, len(gal)):
    for h in range(0,1):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] #/ (1. + zem[h])
        flux = data['norm']
        fx = data['fx']
        var = data['var']

        ## Error Bar Calculations

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * np.sqrt(var[i])/fx[i]

                                                  ##  What connects our code to linetools code ##
                                                  
        xspec=XSpectrum1D.from_tuple((wave,flux,sigma))
        

        #Background info
        lt_path = imp.find_module('linetools')[1]
        #xspec = XSpectrum1D.from_file(lt_path+'/spectra/tests/files/UM184_nF.fits')   this was linetools code that we needed to understand in order to get order data to be inserted


        ##Generate Lyman -alpha Lines
        
        mgiitrans = ['MgII 2796', 'MgII 2803']
        #for loop manipulating background info
        abslines = []
        abslines.attrib['z'] = []
        for trans in mgiitrans:
            print(trans, zem[h])
            iline = AbsLine(trans, z=zem[h], vlim=[-3000.,500.]*u.km/u.s)
            #clear_CACHE_LLIST=True
            #iline.z = zem[h]
            #iline.attrib['z'] = zem[h]
            iline.analy['vlim'] = [-3000.,500.]*u.km/u.s
            iline.analy['spec'] = xspec
            # iline.analy['spec'] = xspec  this was for the default code pulled from linetools
            abslines.append(iline)
            iline.attrib['z'] = zem[h]
            # abslines.attrib['z'].append(iline.z = zem[h])
            abslines.attrib['z'].append(iline.attrib['z'])
        print('check out my abslines:', abslines)
        
        # abslines.attrib['N'] = 10**14./u.cm**2  # log N
        # abslines.attrib['b'] = 25.*u.km/u.s
        # abslines.attrib['z'] = 2.0

        

        #Generates Voigt Profile
        
        # abslines.analy['spec'] = lsio.readspec('../../linetools/spectra/tests/files/UM184_nF.fits') <-- this is gonna be a problem

        vmodel = abslines.generate_voigt()

        plt_line(vmodel)


        
        #Use self-generated Wavelength  ** This could be the key**

        # N deals with column density
        # b is the Doppler width of the voigt profile
        
        abslines.attrib['N'] = 10**17.5/u.cm**2   # ** Just need to figure out both these parameters apparently
        abslines.attrib['b'] = 300.*u.km/u.s

        # wave = np.linspace(3644, 3650,100)*u.AA   # ** insert our velocity ranges here **
        wave = np.linspace(500, 2500,100)*u.AA
        vmodel2 = abslines.generate_voigt(wave=wave)
        plt_line(vmodel2)


        # #Multiple Lines

        # abslin2 = AbsLine('DI 1215')
        # abslin2.attrib['N'] = 10**13./u.cm**2  # log N
        # abslin2.attrib['b'] = 15.*u.km/u.s
        # abslin2.attrib['z'] = 2.0
        # vmodel3 = lav.voigt_from_abslines(wave,[abslin,abslin2])
        # plt_line(vmodel3)

        # #Tau Profiles
        # tau = lav.voigt_from_abslines(wave,abslin,ret='tau')

        # plt.clf()
        # plt.figure(dpi=1200)
        # plt.plot(wave,tau, 'k-', drawstyle='steps-mid', lw=1.5)
        # plt.xlim(3644., 3650.)
        # plt.ylabel('Optical Depth', fontsize=20.)
        # plt.xlabel('Wavelength', fontsize=20.)
        # plt.show()
        # plt.close()

        
        # #Fit the Voigt
        # from astropy.modeling import fitting
        # fitvoigt = lav.single_voigt_model(logN=np.log10(abslines.attrib['N'].value),b=abslines.attrib['b'].value,
        #                         z=abslines.attrib['z'], wrest=abslines.wrest.value, 
        #                         gamma=abslines.data['gamma'].value, f=abslines.data['f'],
        #                          fwhm=3.)
        # fitter = fitting.LevMarLSQFitter()
        # p = fitter(fitvoigt,wave.value,vmodel2.flux.value,weights=1./(np.ones(len(wave))*0.1))
        # print(p)

        # plt.clf()
        # plt.plot(wave, vmodel2.flux, 'k-',drawstyle='steps')
        # plt.plot(wave,p(wave.value), 'g-')
        # #xdb.xplot(wave,model.flux,p(wave.value))
        # plt.show()
        # plt.close()
        # # note the bad pixels at either end of the spectrum - this is a bug! It can't be 
        # # worked around for now by making sure that the spectrum is much larger than the region
        # # of interest, and discarding the first and last ~10 pixels.




        
