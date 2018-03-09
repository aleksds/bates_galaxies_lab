# Kwamae Delva
# March 9, 2018
#Modifies Aleks original code with our data
# https://github.com/linetools/linetools/blob/master/docs/AbsLine_examples.rst --> doesn't work
# https://github.com/linetools/linetools/blob/master/docs/examples/AbsComponent_ColumnDensities.ipynb --> what the below is based on

#import sys
#sys.path.append('/Users/aleks/github/linetools/')

# suppress warnings for these examples
import warnings
warnings.filterwarnings('ignore')

# import
import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools import spectralline as ltsp
from linetools.spectra.xspectrum1d import XSpectrum1D

import numpy as np
try:
    import seaborn as sns; sns.set(context="notebook",font_scale=2)
except:
    pass

from pkg_resources import resource_filename
from scipy import integrate

from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.lists.linelist import LineList

ism = LineList('ISM')

# read in the spectrum
# xspec = XSpectrum1D.from_file(resource_filename('linetools','/spectra/tests/files/UM184_nF.fits'))


 ##  What connects our code to linetools code ##
                                                  
xspec=XSpectrum1D.from_tuple((wave,flux,sigma))


        

# generate AbsLines
Mgtrans = ['MgII2796', 'MgII 2803', 'MgI 2852']
FeIItrans = ['Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']

abslines_Mg = []
for trans in Mgtrans:
    iline = AbsLine(trans,z=zem[h], linelist=ism)
    iline.limits.set([-3000.,500.]*u.km/u.s) # vlim
    iline.analy['spec'] = xspec
    abslines.append(iline)
#
print(abslines_Mg)

abslines_Fe = []
for trans in FeIItrans:
    iline = AbsLine(trans,z=zem[h], linelist=ism)
    iline.limits.set([-3000.,500.]*u.km/u.s) # vlim
    iline.analy['spec'] = xspec
    abslines_Fe.append(iline)
#
print(abslines_Fe)

# Generate the Component
abscomp_Mg = lt_abscomp.AbsComponent.from_abslines(abslines_Mg)
abscomp_Fe = lt_abscomp.AbsComponent.from_abslines(abslines_Fe)

try:
    sns.set(context="notebook",font_scale=2)
except:
    pass
abscomp_Mg.stack_plot()
abscomp_Fe.stack_plot()

# It sucks that stack_plot is the only plotting function for this type of code, but it does allow for direct comparison within a galaxy

# Synthesize / Measure AODM Column Densities
abscomp.synthesize_colm(redo_aodm=True)

print(abscomp_Mg.logN)
print(abscomp_Mg.sig_logN)
print(abscomp_Fe.logN)
print(abscomp_Fe.sig_logN)

for iline in abscomp._abslines:
    print(iline.wrest, iline.attrib['flag_N'], iline.attrib['logN'], iline.attrib['sig_logN'])

# Apparent Column Density Plot
abscomp.plot_Na()



# This code shows relationship between ????????

# Curve of Growth
def ftau_intgrnd(x,tau0=0.1):
    return 1 - np.exp(-tau0 * np.exp(-x**2))

neval = 10000
lgt = np.linspace(-3, 9, neval)
all_tau0 = 10.**lgt
Ftau = np.zeros(neval)
for jj,tau0 in enumerate(all_tau0):
    Ftau[jj], ferr = integrate.quad(ftau_intgrnd, 0, np.inf, args=(tau0,))

# Damped limit (not accurate enough)
damp_lgt = np.linspace(6, 10, 100)
damp_tau0 = 10.**damp_lgt
damp_Ftau = np.sqrt(np.log(damp_tau0))

import matplotlib.pyplot as plt

plt.plot(lgt, Ftau, damp_lgt, 1.015*damp_Ftau)
plt.show()

# Perform and Plot
abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)
COG_dict = abscomp.cog(redo_EW=True, show_plot=True)

plt.show()

print(COG_dict)


