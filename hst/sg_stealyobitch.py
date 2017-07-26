# sophia gottlieb
# 20170711
# modified from cl_std.py by charlie lipscomb
import numpy as np
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg
from mmap import mmap,ACCESS_READ
from xlrd import open_workbook

# defines a galaxy to have a name, redshift z, and initial x and y coordinates
class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        

res = ['fine','coarse']
wb = open_workbook(res[1]+'data.xlsx')

for sheet in wb.sheets():
    numr = sheet.nrows
    Galaxies = []

    for row in range(1,numr):
        Galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))
        

#define a function to plot "postage stamp" images
def plot_image(posx,posy, count, prevstds):
    stamp = data[int(round(posy-dy)):int(round(posy+dy)), int(round(posx-dx)):int(round(posx+dx))]

    count = count +1

    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    xavg = np.average([x1,x2,x3])
    yavg = np.average([y1,y2,y3])
    xstd = np.std([x1,x2,x3])
    ystd = np.std([y1,y2,y3])

    # RECURSION BITCH
    if count < 100 and (xstd > 0.1 or ystd > 0.1) and (xstd <= prevstds[0] and ystd <= prevstds[1]):
        print(count, posx-dx+xavg, posy-dy+yavg, xstd, ystd)
        return plot_image(posx-dx+xavg, posy-dy+yavg, count, [xstd,ystd])
    else:
        return posx-dx+xavg, posy-dy+yavg, 1/(xstd**2), 1/(ystd**2), count


# define the directory that contains the images
dir = os.environ['HSTDIR']


dx = 5
dy = 5

# create a PDF file for the plots
meth = ['centroid_com','centroid_1dg','centroid_2dg']
foureightone = [4,8,1]
filters = ['F475W','F814W','F160W']
colors = ['blue','green','red']# or these same things with lengal and len fil
xcen = np.zeros([len(Galaxies)])
ycen = np.zeros([len(Galaxies)])

xstds = np.zeros([len(Galaxies)])
ystds = np.zeros([len(Galaxies)])

for j in range(0,1):
#for j in range(0,len(Galaxies)):
    #xstds[j] = 1
    #ystds[j] = 1
        mxs = [0,0,0]
        mys = [0,0,0]
        mstdx = [0,0,0]
        mstdy = [0,0,0]
        for i in range(0,len(filters)):
            file = glob.glob(dir+Galaxies[j].name+'/'+res[1]+'/'+str(filters[i])+'/final_'+str(filters[i])+'*sci.fits')
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header

            mxs[i], mys[i], mstdx[i], mstdy[i], count = plot_image(Galaxies[j].x,Galaxies[j].y, 0, [1000, 1000])
            
        xcen[j] = np.average(mxs, weights = mstdx)
        ycen[j] = np.average(mys, weights = mstdy)
        

