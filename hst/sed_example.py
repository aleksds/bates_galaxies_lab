# Aleks Diamond-Stanic
# 20190506
# the goal is to plot an example spectral energy distribution

# import relevant modules
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Here is photometric information for one galaxy.
# This includes the following:
# [0] flux in units of maggies
# [1] inverse variance in units of maggies
# [2] effective wavelength in units of microns

phot = dict(
        FUV=(1.82491e-09, 1.15115e+19, 0.1528),
        NUV=(8.06441e-09, 2.25963e+19, 0.2271),
        u=(1.27666e-08, 1.53319e+18, 0.3543),
        g=(1.59991e-08, 7.47418e+18, 0.4770),
        r=(3.15573e-08, 2.18797e+18, 0.6231),
        i=(3.63049e-08, 1.53877e+18, 0.7625),
        z=(4.14564e-08, 2.71207e+17, 0.9134),
        W1=(9.40651e-08, 5.61366e+17, 3.368),
        W2=(1.02882e-07, 1.38784e+17, 4.618),
        W3=(5.44324e-07, 8.12757e+14, 12.082),
        W4=(1.38524e-06, 1.90498e+13, 22.194))

# use the information above to define wavelength and flux arrays
wave = np.array([phot[filt][2] for filt in phot.keys()])
flux = np.array([phot[filt][0] for filt in phot.keys()])

# made an output PDF file
file = 'sed_example.pdf'

with PdfPages(file) as pdf:

    # create a figure using matplotlib
    fig = plt.figure()

    # plot wavelength on the x axis and flux on the y axis
    plt.scatter(wave, flux)

    # set parameters of the y axis
    plt.yscale('log')
    plt.ylim([1e-10,1e-5])
    plt.ylabel('Flux [maggies]')

    # set parameters of the x axis
    plt.xscale('log')
    plt.xlim([0.1, 30.])
    plt.xlabel('Wavelength [microns]')

    # save the figure
    pdf.savefig()
    plt.close()

# open the pdf file
os.system('open %s &' % file)
