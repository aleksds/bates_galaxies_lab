# Aleks Diamond-Stanic, 20170623
# goal is to test out plotting Voigt profiles using Voigt1D from astropy

import numpy as np
import os
from astropy.modeling.models import Voigt1D
import matplotlib.pyplot as plt
from astropy.modeling import fitting
from matplotlib.backends.backend_pdf import PdfPages

filename = 'ad_voigt_test.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    x = np.arange(-3500, 500, 1)
    v1 = Voigt1D(x_0=-1000, amplitude_L=-1, fwhm_L=100, fwhm_G=100)

    ax.plot(x, v1(x)+1.)
    plt.ylabel('Continuum-normalized Flux')
    plt.xlabel('Velocity [km/s]')

    pdf.savefig()
    plt.close()
            
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    voi = Voigt1D(amplitude_L=-0.5, x_0=0.0, fwhm_L=5.0, fwhm_G=5.0)
    xarr = np.linspace(-5.0,5.0,num=40)
    yarr = voi(xarr)
    
    voi_init = Voigt1D(amplitude_L=-1.0, x_0=1.0, fwhm_L=5.0, fwhm_G=5.0)
    fitter = fitting.LevMarLSQFitter()
    voi_fit = fitter(voi_init, xarr, yarr)
    print(voi_fit)

    ax.plot(xarr,yarr)
    ax.plot(xarr,voi_fit(xarr))

    pdf.savefig()
    plt.close()

os.system("evince %s &" % filename)    
