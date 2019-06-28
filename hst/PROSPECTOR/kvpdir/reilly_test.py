#testing reilly's code

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


fraction = np.array([7 / 60, 3 / 59, 8 / 59, 3 / 59, 9 / 59, 8 / 59, 9 / 58, 16 / 59, 23 / 60])
unc = np.array([np.sqrt(7) / 60, np.sqrt(3) / 59, np.sqrt(8) / 59, np.sqrt(3) / 59, np.sqrt(9) / 59, np.sqrt(8) / 59,
                np.sqrt(9) / 58,
                np.sqrt(16) / 59, np.sqrt(23) / 60])
sigma_sfr = np.array([0.005, 0.001, 0.0015, 0.002, 0.003, 0.0036, 0.0047, 0.007, 0.01])

filename = 'reilly_test.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.scatter(sigma_sfr, fraction)
    # plt.errorbar(sigma_sfr, fraction, yerr=unc)
    plt.xlabel('SFR Surface Density')
    plt.ylabel('Fraction')
    #plt.xlim(-5, -1)
    #plt.ylim(0, 1)

    pdf.savefig()
    plt.close

os.system('open %s &' % filename)