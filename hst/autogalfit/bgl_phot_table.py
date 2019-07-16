# Aleks Diamond-Stanic 2019 July 16
# the goal is to read in photometry from a table and then print out colors and uncertainties

import numpy as np
from astropy.io import ascii

data = ascii.read('bgl_phot.dat')

ngal = len(data)

for i in range(0, ngal):
    m475 = data['m475'][i] - data['ebv'][i] * 3.248
    m814 = data['m814'][i] - data['ebv'][i] * 1.536
    m160 = data['m160'][i] - data['ebv'][i] * 0.512
    #uvis_color = data['m475'][i] - data['m814'][i]
    uvis_color = m475 - m814
    uvis_unc = np.sqrt(data['u475'][i]**2 + data['u814'][i]**2)
    #uvir_color = data['m814'][i] - data['m160'][i]
    uvir_color = m814 - m160
    uvir_unc = np.sqrt(data['u814'][i]**2 + data['u160'][i]**2)
    #print(data['Galaxy'][i], "{0:.3f}".format(uvis_color), "{0:.3f}".format(uvis_unc), "{0:.3f}".format(uvir_color), "{0:.3f}".format(uvir_unc))
    print(data['Galaxy'][i]+' & $'+\
              "{0:.3f}".format(m475)+'\pm'+"{0:.3f}".format(data['u475'][i])+'$ & $'+\
              "{0:.3f}".format(m814)+'\pm'+"{0:.3f}".format(data['u814'][i])+'$ & $'+\
              "{0:.3f}".format(m160)+'\pm'+"{0:.3f}".format(data['u160'][i])+'$ & $'+\
              "{0:.3f}".format(uvis_color)+'\pm'+"{0:.3f}".format(uvis_unc)+'$ & $'+\
              "{0:.3f}".format(uvir_color)+'\pm'+"{0:.3f}".format(uvir_unc)+'$ '+'\\'+'\\')

