# Code by Charles Lipscomb with lots of code sampled from Sophia Gottlieb and Aleksandar Diamond-Stanic
# Quick description: This code compares centroid values based on three routines from photutils: centroid_com, centroid_1dg, centroid_2dg.
# Future developments: Could turn this into a function that returns the coordinates of the best centroid location.
# import relevant Python modules
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
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
Galaxies = [['J0826', 3629, 4154], ['J0901', 3934, 4137], ['J0905', 3387, 3503], ['J0944', 3477., 3405.], ['J1107', 3573, 3339.], ['J1219', 3803., 4170.], ['J1341', 3884, 4165], ['J1506', 4147., 3922.], ['J1558', 3787., 4186.], ['J1613', 4175., 3827.], ['J2116', 3567, 3436], ['J2140', 4067, 4054]]
dx = 5
dy = 5
methods = ['centroid_com','centroid_1dg','centroid_2dg']
foureightone = [4,8,1]
filters = ['F475W','F814W','F160W']
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
bestxs = np.zeros(12)
bestys = np.zeros(12)

#calculation of best value for centroids
for w in range(0,len(Galaxies)):
    mxs = [0,0,0]
    mys = [0,0,0]
    for i in range(0, len(filters)):
        file = glob.glob(dir+galaxies[w]+'*/coarse/'+filters[i]+'/final*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']
        gain[i] = header[i]['CCDGAIN']
        RN[i] = header[i]['READNSEA']
        mxs[i], mys[i] = plot_image_monster(Galaxies[w][1], Galaxies[w][2])
    bestxs[w] = np.average(mxs)
    bestys[w] = np.average(mys)

#make pdf for plots and plot
filename = 'cl_centroid.pdf'
with PdfPages(filename) as pdf:
    for j in range(0,len(Galaxies)):
        fig = plt.figure()
        plt.text(Galaxies[j][1]+0.36,Galaxies[j][2]-.05, xstds[j])
        plt.text(Galaxies[j][1]+0.2,Galaxies[j][2]-.05, 'sigma_x =')
        plt.text(Galaxies[j][1]+0.36,Galaxies[j][2]-0.1, ystds[j])
        plt.text(Galaxies[j][1]+0.2,Galaxies[j][2]-.1, 'sigma_y =')
        
        plt.scatter(Galaxies[j][1],Galaxies[j][2], label='brightness center',color='black')
        plt.scatter(bestxs[j],bestys[j], label='BEST',color='gold')
        for i in range(0,len(filters)):
            file = glob.glob(dir+galaxies[j]+'*/coarse/'+filters[i]+'/final*sci.fits')
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
            stamp = data[round(Galaxies[j][2]-dy):round(Galaxies[j][2]+dy), round(Galaxies[j][1]-dx):round(Galaxies[j][1]+dx)]
            coor = plot_image()
            for m in range(0,len(methods)):
                if i==0:
                    plt.scatter(coor[0][m]+Galaxies[j][1]-dx,coor[1][m]+Galaxies[j][2]-dy,label = methods[m],color=colors[m])
                else:
                    plt.scatter(coor[0][m]+Galaxies[j][1]-dx,coor[1][m]+Galaxies[j][2]-dy,color=colors[m],marker=markers[i])
            plt.xlabel('x')
            plt.ylabel('y')
            plt.xlim(Galaxies[j][1]-.5,Galaxies[j][1]+.5)
            plt.ylim(Galaxies[j][2]-.5,Galaxies[j][2]+.5)
            plt.title(Galaxies[j][0] + ' Centroid Calculations')
            legend = plt.legend(loc='upper left')
            plt.text(Galaxies[j][1]+0.2,Galaxies[j][2]+0.15, r'Circle = F475W')
            plt.text(Galaxies[j][1]+0.2,Galaxies[j][2]+0.1, r'Plus = F814W')
            plt.text(Galaxies[j][1]+0.2,Galaxies[j][2]+0.05, r'Diamond = F160W')
            
        pdf.savefig()
        plt.close()
    os.system('open %s &' % filename)

