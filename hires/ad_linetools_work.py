# Aleks Diamond-Stanic
# Feb 16, 2018
#
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
xspec = XSpectrum1D.from_file(resource_filename('linetools','/spectra/tests/files/UM184_nF.fits'))

# generate AbsLines
SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
abslines = []
for trans in SiIItrans:
    iline = AbsLine(trans,z=2.92939, linelist=ism)
    iline.limits.set([-250.,80.]*u.km/u.s) # vlim
    iline.analy['spec'] = xspec
    abslines.append(iline)
#
print(abslines)

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


