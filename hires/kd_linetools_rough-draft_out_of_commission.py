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

try:
    import seaborn as sns; sns.set(context="notebook",font_scale=2)
except:
    pass

warnings.filterwarnings('ignore')

#Background info
lt_path = imp.find_module('linetools')[1]
xspec = XSpectrum1D.from_file(lt_path+'/spectra/tests/files/UM184_nF.fits')
names = ['Mg II 2796', 'Mg II 2803', 'Mg I 2852']

#for loop manipulating background info
abslines = []
for trans in names:
    iline = AbsLine(trans)
    clear_CACHE_LLIST=True
    iline.attrib['z'] = .603
    iline.analy['vlim'] = [-3000.,500.]*u.km/u.s
    iline.analy['spec'] = xspec
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
