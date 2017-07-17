#so_finedata.py
# the goal of this code is to access the fine data for the galaxies that have fine data

import numpy as np
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg
from mmap import mmap,ACCESS_READ
from xlrd import open_workbook

dir = '/Volumes/physics/linux-lab/data/hst/'


# defines a galaxy to have a name, redshift z, and initial x and y coordinates
class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        
wb = open_workbook('finedata.xlsx')

for sheet in wb.sheets():
    numr = sheet.nrows
    galaxies = []

    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))
wavelengths = ['F475W','F814W']
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
for w in range(0,len(galaxies)):
    print(galaxies[w].name)        
    for i in range(0, len(wavelengths)):
        file = glob.glob(dir+galaxies[w].name+'/fine/'+wavelengths[i]+'/final*sci.fits')
            
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
