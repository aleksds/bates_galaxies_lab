import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from scipy import integrate
from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM 
from mmap import mmap,ACCESS_READ
from xlrd import open_workbook
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

#function that runs all 3 methods
def plot_image():
    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    return np.array([[x1,x2,x3] , [y1,y2,y3]])


# define the directory that contains the images
dir = os.environ['HSTDIR']

# parameters for the initial guess at the galaxy centroid and the size
# pix coordinates on entire ds9 image
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J2116', 'J2140']
xcen = [7260,7870,6775,6956,7148,7609,7769,8296,7576,7134,8134]
ycen = [8310,8276,7009,6813,6680,8343,8332,7846,8374,6872,8110]

dx = 5
dy = 5
methods = ['centroid_com','centroid_1dg','centroid_2dg']
foureightone = [4,8,1]
filters = ['F475W','F814W']
colors = ['blue','green','red']
markers = ('o','+','d')

#results from cl_std.py
xstds = [0.06,  0.03,  0.09,  0.06,  0.07, 0.09,  0.08,  0.04,  0.04,  0.16, 0.06,  0.06]
ystds = [0.07,  0.05,  0.11 ,  0.08,  0.08, 0.10,  0.06,  0.08,  0.08,  0.11, 0.14,  0.10]

#define function that doesnt use center of mass method
def plot_image_monster(posx,posy):
    stamp = data[i][int(round(posy-dy)):int(round(posy+dy)), int(round(posx-dx)):int(round(posx+dx))]
    std = np.std(stamp[stamp==stamp])
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    xavg = ((x2-dx+posx)+(x3-dx+posx))/2
    yavg = ((y2-dy+posy)+(y3-dy+posy))/2
    xstd = np.std([x2,x3])
    ystd = np.std([y2,y3])
    return xavg, yavg

#empty lists to be refreshed and filled
data = [0 for x in range(len(filters))]
header = [0 for x in range(len(filters))]
fnu = [0 for x in range(len(filters))]
exp = [0 for x in range(len(filters))]
gain = [0 for x in range(len(filters))]
dark = [0.0022,0.0022,0.045]
RN = [0 for x in range(len(filters))]
bestxs = np.zeros(11)
bestys = np.zeros(11)

#calculation of best value for centroids
for w in range(0,len(galaxies)):
    mxs = [0,0]
    mys = [0,0]
    for i in range(0, len(filters)):
        file = glob.glob(dir+galaxies[w]+'*/fine/'+filters[i]+'/final*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']
        gain[i] = header[i]['CCDGAIN']
        RN[i] = header[i]['READNSEA']
        mxs[i], mys[i] = plot_image_monster(xcen[w], ycen[w])
    bestxs[w] = np.average(mxs)
    bestys[w] = np.average(mys)
