import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from os.path import isdir, join
import csv
import math
import marvin


filename = 'map'
with PdfPages(filename) as pdf:
    my_cube = marvin.tools.Maps('7443-12703')
    ha = my_cube.emline_gflux_ha_6564
    fig, ax = ha.plot()

    pdf.savefig()
    plt.close()

os.system("open %s &" % filename)