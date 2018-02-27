# Kwamae Delva
#Code using linetools to make a better voigt profile plot
#Code implementing column density





from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.lists.linelist import LineList


import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import astropy.units as u
from linetools import spectralline as ltsp
from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.analysis import voigt as lav
import imp
import warnings
from astropy.io import ascii
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting

from pkg_resources import resource_filename
from scipy import integrate

try:
    import seaborn as sns; sns.set(context="notebook",font_scale=2)
except:
    pass

warnings.filterwarnings('ignore')

ism = LineList('ISM')

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
    for h in range(1,2):
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


        mgiitrans = ['MgII 2796', 'MgII 2803']
        #for loop manipulating background info
        abslines = []
        for trans in mgiitrans:
            print(trans, zem[h])
            iline = AbsLine(trans, z=zem[h], linelist=ism)
            iline.limits.set([-3000.,400.]*u.km/u.s) # vlim
            #clear_CACHE_LLIST=True
            #iline.z = zem[h]
            #iline.attrib['z'] = zem[h]
            #iline.analy['vlim'] = [-3000.,500.]*u.km/u.s
            iline.analy['spec'] = xspec
            # iline.analy['spec'] = xspec  this was for the default code pulled from linetools
            abslines.append(iline)
        print('check out my abslines:', abslines)

            #Wavelengths of interest
            #MgII2796wrest=2796.3542699
            #MgII2803wrest=2803.5314853
            #MgI2852wrest=2852.96328

            #AbsLine=[MgII2796wrest, MgII2803wrest, MgI2852wrest]


        abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)
        #try:
        #    sns.set(context="notebook",font_scale=2)
        #except:
        #    pass
        #Generates components for column density plot
        abscomp.stack_plot(vlim=[-3000.,500.]*u.km/u.s)

        #Produces numbers for plot generation for each wavelength
        zlim=[zem[h]-0.01, zem[h]+0.002]
        kwargs={}
        kwargs['zlim'] = zlim
        abscomp.synthesize_colm(redo_aodm=True)#, **kwargs)#limits=[-3000.,500.]*u.km/u.s)#, vlim=[-3000.,500.]*u.km/u.s)
        print('give me some column density:', abscomp.logN)
        for iline in abscomp._abslines:
            print(iline.wrest, iline.attrib['flag_N'], iline.attrib['logN'], iline.attrib['sig_logN'])

            #Plots Apparent Column Density
            abscomp.plot_Na()

        # adding in code from ad_linetools_exvo.py (clearly need multiple components)
        fitvoigt = lav.single_voigt_model(logN=np.log10(abscomp.attrib['N'].value),b=100,#b=abscomp.attrib['b'].value,
                                z=abslines[0].z, wrest=abslines[0].wrest.value, 
                                gamma=abslines[0].data['gamma'].value, f=abslines[0].data['f'],
                                 fwhm=3.)

        fitter = fitting.LevMarLSQFitter()


        p = fitter(fitvoigt,xspec.wavelength.value,xspec.flux.value,weights=1./(np.ones(len(xspec.wavelength.value))*0.1))
        print(p)

        plt.plot(xspec.wavelength.value, xspec.flux.value, 'k-',drawstyle='steps')
        plt.plot(xspec.wavelength.value,p(xspec.wavelength.value), 'g-')
        plt.xlim(4050., 4100.)
        plt.ylim(0., 1.1)

        plt.show()
        plt.close()
