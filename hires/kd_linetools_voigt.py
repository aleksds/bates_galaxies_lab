# Kwamae Delva
#Code using linetools to make a better voigt profile plot
#Code implementing column density


import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from scipy import integrate
import astropy.units as u
from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectralline import AbsLine
from linetools.spectra.xspectrum1d import XSpectrum1D
import imp
import warnings
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting

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
## What does this do?????????????
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

                                                  ##  What connects our code to linetools code ##
                                                  
        xspec=XSpectrum1D.from_tuple((wave,flux,sigma))
        

        #Background info
        lt_path = imp.find_module('linetools')[1]
        #xspec = XSpectrum1D.from_file(lt_path+'/spectra/tests/files/UM184_nF.fits')   this was linetools code that we needed to understand in order to get order data to be inserted

        #for loop manipulating background info
        abslines = []
        for trans in names:
            iline = AbsLine(trans)
            clear_CACHE_LLIST=True
            iline.attrib['z'] = zem[h]
            iline.analy['vlim'] = [-3000.,500.] #*u.km/u.s
            iline.analy['spec'] = xspec
            # iline.analy['spec'] = xspec  this was for the default code pulled from linetools
            abslines.append(iline)




            #Wavelengths of interest
            MgII2796wrest=2796.3542699
            MgII2803wrest=2803.5314853
            MgI2852wrest=2852.96328

            AbsLine=[MgII2796wrest, MgII2803wrest, MgI2852wrest]


            abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)
            try:
                sns.set(context="notebook",font_scale=2)
            except:
                pass
            #Generates components for column density plot
            abscomp.stack_plot()

            #Produces numbers for plot generation for each wavelength
            abscomp.synthesize_colm(redo_aodm=True)
            abscomp.logN
            for iline in abscomp._abslines:
                print(iline.wrest, iline.attrib['flag_N'], iline.attrib['logN'], iline.attrib['sig_logN'])

                #Plots Apparent Column Density
                abscomp.plot_Na()
