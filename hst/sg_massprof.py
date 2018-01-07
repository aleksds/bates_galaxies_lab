#sophia gottlieb
#dec 29 2017
# mass profiles

import os
cwd = os.getcwd()
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
import pandas as pd

### THINGS FOR READING IN GALAXY DATA FROM galaxydata.xlsx

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

from mmap import mmap,ACCESS_READ
from xlrd import open_workbook

# for writing into the xl file
import xlwt

#Define Galaxy class to hold name, redshift, x and y positions, and
# Have functions for luminosity distance in cms and rad to kpc conversion
# Okay so i checked and the functions don't work and I hate them. Anyway,
# I saved them as fields and it's working fine now
class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = 0.05/(conv.arcsec_per_kpc_proper(self.z)/u.arcsec*u.kpc)

radii = np.arange(40)+1
area = [0 for x in range(len(radii))]

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

        
# Grabbing and filling galaxy data
res = ['fine','coarse']
t=0

type = ['convolved_image', 'final']

#Grab galaxy data from excel
wbo = open_workbook('galaxydata.xlsx')

for sheet in wbo.sheets():
    numr = sheet.nrows
    galaxies = []
    
    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))
count = -1

pltcolors = ['b','g','r','c','m','y','k']
pltmarker = ['o','*','v','P','h','d']

datadir = '/Users/sgottlie/Downloads/subset_20171220/'
with PdfPages('sg_massprofiles_UVIS.pdf') as pdf:
    for root, dirs, filenames in os.walk(datadir):
        for f in filenames:
            if "UVIS" in f and ".xls" in f:
                count += 1
                gal = galaxies[count]
                data = pd.read_excel(datadir+f)
                print(datadir+f)

                flux = pd.read_csv(cwd + '/PROSPECTOR/Photometry/coarse_final/'+gal.name[:5]+'_uvis.txt', sep = "\t", header = 0)

                light = flux['Unnamed: 3']
                totallight = np.sum(light)
                totalmass = np.sum(data['best mass'])
                mlrscalar = totalmass / totallight

                
                kpcArea = area*gal.radToKpc**2
                fig = plt.figure(figsize = (10,8))
                
                ax = fig.add_subplot(2,1,1)
                plt.title('Mass vs Radii for ' + gal.name)
                #plt.xlabel('Radii (kpc)')
                plt.ylabel('Mass (solar masses)')
                
                plt.semilogy(radii*gal.radToKpc, light*mlrscalar,linestyle = 'None',  marker = 'o', color = 'g', label = 'Light Profile [F814W]')
                plt.semilogy(radii*gal.radToKpc, data['best mass'],linestyle = 'None',  marker = '*', color = 'b', label ='Stellar Mass Profile')
                plt.legend()


                ax = fig.add_subplot(2,1,2)
                plt.title('Mass Density vs Radii for '+ gal.name)
                plt.xlabel('Radii (kpc)')
                plt.ylabel('Mass Density (solar masses / kpc^2)')
                
                plt.semilogy(radii*gal.radToKpc, (light*mlrscalar)/kpcArea,linestyle = 'None',  marker = 'o', color = 'g', label = 'Surface Brightness Profile [F814W]')
                plt.semilogy(radii*gal.radToKpc, data['best mass']/kpcArea,linestyle = 'None',  marker = '*', color = 'b', label = 'Stellar Mass Surface Density Profile')
                plt.legend()
                
                #plt.close()

                
                pdf.savefig()
                plt.close()

