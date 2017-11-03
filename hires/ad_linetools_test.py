# imports
try:
    import seaborn as sns; sns.set(context="notebook",font_scale=2)
except:
    pass

from scipy import integrate
import astropy.units as u

from linetools.isgm import abscomponent as lt_abscomp
from linetools.spectralline import AbsLine
from linetools.spectra.xspectrum1d import XSpectrum1D
#
import imp
lt_path = imp.find_module('linetools')[1]


MgIItrans = ['MgII 2796', 'MgII 2803']

abslines = []
for trans in MgIItrans:
    iline = AbsLine(trans)
    iline.attrib['z'] = 2.92939
    iline.analy['vlim'] = [-250.,80.]*u.km/u.s
    #iline.analy['spec'] = xspec
    abslines.append(iline)
#
abslines

abscomp = lt_abscomp.AbsComponent.from_abslines(abslines)
