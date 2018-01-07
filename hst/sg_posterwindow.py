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
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM 
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

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
dir = '/Users/sgottlie/Desktop/linux-lab/'


# define a function to plot "postage stamp" images
def plot_image():
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='Greys')
    plt.tick_params(axis='both', which='major', labelsize=8)
    
filters = ['F475W','F814W','F160W']
res = ['fine','coarse']
type = ['convolved_image', 'final']


alldata = []

with PdfPages('poster.pdf') as pdf:
    for g in range(0,1):
        #for g, gal in enumerate(galaxies):
        gal = galaxies[g]
        fig = plt.figure()
        plt.suptitle(gal.name)
        for f, fil in enumerate(filters):
            file = glob.glob(dir+gal.name+'*/'+res[1]+'/'+fil+'/'+type[1]+'_'+fil+'*sci.fits')

            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
            width = 100
            dy = width
            dx = width

            stamp = data[round(gal.y-dy):round(gal.y+dy), round(gal.x-dx):round(gal.x+dx)]
            #ax = fig.add_subplot(len(galaxies),len(filters), 1+g*4+f)
            ax = fig.add_subplot(1,len(filters)+1, 1+f)
            plt.axis('off')
            plot_image()
            plt.title(fil)
            alldata.append(stamp)


        alldata = np.reshape(alldata, (200,200,3))
        alldata = alldata - np.min(alldata)
        alldtat = alldata/np.max(alldata)
        bx = fig.add_subplot(1, len(filters)+1, 4)
        plt.imshow(alldata)    
        plt.axis('off')
        pdf.savefig()
        plt.close()
                             
