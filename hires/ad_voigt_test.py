# Aleks Diamond-Stanic, 20170623
# goal is to test out plotting Voigt profiles using Voigt1D from astropy

import numpy as np
import os
from astropy.modeling.models import Voigt1D
import matplotlib.pyplot as plt
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
            
os.system("evince %s &" % filename)    
