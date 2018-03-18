# Aleks Diamond-Stanic
# Feb 16, 2018
#
# https://github.com/linetools/linetools/blob/master/docs/AbsLine_examples.rst --> doesn't work
# https://github.com/linetools/linetools/blob/master/docs/examples/AbsComponent_ColumnDensities.ipynb --> what the below is based on

# Kwamae Delva
#Code using linetools to make a better voigt profile plot
#Code implementing column density

import sys
sys.path.append('/Users/kdelva/github/linetools/')



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

        FeIItrans = ['Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']
        abslines = []
        for trans in FeIItrans:
            iline = AbsLine(trans,z=zem[h], linelist=ism)
            iline.limits.set([-3000.,400.]*u.km/u.s) # vlim
            iline.analy['spec'] = xspec
            abslines.append(iline)

        # Generate the Component
        abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)

        try:
            sns.set(context="notebook",font_scale=2)
        except:
            pass
        abscomp.stack_plot()

        # Synthesize / Measure AODM Column Densities
        abscomp.synthesize_colm(redo_aodm=True)

        print(abscomp.logN)
        print(abscomp.sig_logN)

        for iline in abscomp._abslines:
            print(iline.wrest, iline.attrib['flag_N'], iline.attrib['logN'], iline.attrib['sig_logN'])

        # Apparent Column Density Plot
        abscomp.plot_Na()




# import sys
# sys.path.append('/Users/kdelva/github/linetools/')

# # suppress warnings for these examples
# import warnings
# warnings.filterwarnings('ignore')

# # import
# import astropy.units as u
# from linetools.spectralline import AbsLine, SpectralLine
# from linetools import spectralline as ltsp
# from linetools.spectra.xspectrum1d import XSpectrum1D

# import numpy as np
# try:
#     import seaborn as sns; sns.set(context="notebook",font_scale=2)
# except:
#     pass

# from pkg_resources import resource_filename
# from scipy import integrate

# from linetools.isgm import abscomponent as lt_abscomp
# from linetools.spectra.xspectrum1d import XSpectrum1D
# from linetools.lists.linelist import LineList

# ism = LineList('ISM')

# # read in the spectrum
# xspec = XSpectrum1D.from_file(resource_filename('linetools','/spectra/tests/files/UM184_nF.fits'))

# # # generate AbsLines
# # SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
# # abslines = []
# # for trans in SiIItrans:
# #     iline = AbsLine(trans,z=2.92939, linelist=ism)
# #     iline.limits.set([-250.,80.]*u.km/u.s) # vlim
# #     iline.analy['spec'] = xspec
# #     abslines.append(iline)
# #

# FeIItrans = ['Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']
# abslines = []
# for trans in FeIItrans:
#     iline = AbsLine(trans,z=.45900, linelist=ism)
#     iline.limits.set([-3000.,100.]*u.km/u.s) # vlim
#     iline.analy['spec'] = xspec
#     abslines.append(iline)
# print(abslines)

# Generate the Component
abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)

try:
    sns.set(context="notebook",font_scale=2)
except:
    pass
abscomp.stack_plot()

# Synthesize / Measure AODM Column Densities
abscomp.synthesize_colm(redo_aodm=True)

print(abscomp.logN)
print(abscomp.sig_logN)

for iline in abscomp._abslines:
    print(iline.wrest, iline.attrib['flag_N'], iline.attrib['logN'], iline.attrib['sig_logN'])

# Apparent Column Density Plot
abscomp.plot_Na()

# # # Curve of Growth
# # def ftau_intgrnd(x,tau0=0.1):
# #     return 1 - np.exp(-tau0 * np.exp(-x**2))

# # neval = 10000
# # lgt = np.linspace(-3, 9, neval)
# # all_tau0 = 10.**lgt
# # Ftau = np.zeros(neval)
# # for jj,tau0 in enumerate(all_tau0):
# #     Ftau[jj], ferr = integrate.quad(ftau_intgrnd, 0, np.inf, args=(tau0,))

# # # Damped limit (not accurate enough)
# # damp_lgt = np.linspace(6, 10, 100)
# # damp_tau0 = 10.**damp_lgt
# # damp_Ftau = np.sqrt(np.log(damp_tau0))

# # import matplotlib.pyplot as plt

# # plt.plot(lgt, Ftau, damp_lgt, 1.015*damp_Ftau)
# # plt.show()

# # # Perform and Plot
# # abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)
# # COG_dict = abscomp.cog(redo_EW=True, show_plot=True)

# # plt.show()

# # print(COG_dict)


