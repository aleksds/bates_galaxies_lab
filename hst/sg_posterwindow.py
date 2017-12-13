import numpy as np
from astropy.io import fits
import glob
import os
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D

from xlrd import open_workbook

class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc

# Grabbing and filling galaxy data
cf = ['coarse','fine']

# label describing image processing
rc = ['final','convolved_image']
# matching file name ending whoooooooot.
# not confident the h in needed in whooooot. but lets roll with it for now. its a nice h
fex = ['*sci.fits','.fits']

galaxies = []

wbo = open_workbook('galaxydata.xlsx')
for sheet in wbo.sheets():
    numr = sheet.nrows
    
    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))



# define the directory that contains the images
dir = os.environ['HSTDIR']
# define a function to plot "postage stamp" images
def plot_image():
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)


with PdfPAges('poster.pdf') as pdf:
    for w in galaxies:
        fig = plt.figure()
        file = glob.glob(dir+
