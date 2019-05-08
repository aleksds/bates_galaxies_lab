#This is going to be the code for the image comparison between galfit sigma images and our created sigma images

#first recreating our sigma images
#taking code from ad_monster.py

#relevant packages

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
import img_scale

### THINGS FOR READING IN GALAXY DATA FROM galaxydata.xlsx

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

from mmap import mmap, ACCESS_READ
from xlrd import open_workbook

# for writing into the xl file
import xlwt

# Define Galaxy class to hold name, redshift, x and y positions, and
# Have functions for luminosity distance in cms and rad to kpc conversion
# I saved them as fields and it's working fine now
class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        self.lDcm = cosmo.luminosity_distance(self.z) * u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z) * 0.05 / u.arcsec * u.kpc


# setting up final xl file
wb = xlwt.Workbook()
ws = wb.add_sheet('Galaxy Data')

# Grabbing and filling galaxy data
res = ['fine', 'coarse']
t = 0

#type of image we are using
type = ['convolved_image', 'final']

# wbo = open_workbook(res[t]+'data.xlsx')
#loc = ('/Users/adiamond/github/bates_galaxies_lab/hst/galaxydata.xlsx')
loc = ('../finedata.xlsx')
wbo = open_workbook(loc)
# wbo = open_workbook('update.xls')
for sheet in wbo.sheets():
    numr = sheet.nrows
    galaxies = []

    ws.write(0, 0, 'name')
    ws.write(0, 1, 'z')
    ws.write(0, 2, 'x')
    ws.write(0, 3, 'y')

    for row in range(1, numr):
        galaxies.append(Galaxy(sheet.cell(row, 0).value, sheet.cell(row, 1).value, sheet.cell(row, 2).value,
                               sheet.cell(row, 3).value))
        ws.write(row, 0, sheet.cell(row, 0).value)
        ws.write(row, 1, sheet.cell(row, 1).value)

### THINGS FOR CENTROID ###
dx = 5
dy = dx

#plot image function
def plot_image(posx,posy, count, prevstds):
    #We have to redefine the stamp every time, otherwise the program doesn't woek
    stamp = data[i][int(round(posy-dy)):int(round(posy+dy)), int(round(posx-dx)):int(round(posx+dx))]
    #this is just to keep track of how much recursing we do
    count = count +1

    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    xavg = np.average([x2,x3])
    yavg = np.average([y2,y3])
    xstd = np.std([x2,x3])
    ystd = np.std([y2,y3])
    #xavg = np.average([x1,x2,x3])
    #yavg = np.average([y1,y2,y3])
    #xstd = np.std([x1,x2,x3])
    #ystd = np.std([y1,y2,y3])
    #print(count, posx-dx+xavg, posy-dy+yavg, xstd, ystd)
    # RECURSION BITCH limit 100 times, while either std is higher than our 0.1 threshold
    # and as long as the std is getting smaller
    if (xstd + ystd > prevstds[0]+prevstds[1]):
        return posx, posy, prevstds[0]**(-2), prevstds[1]**(-2), count-1
    if count < 100 and (xstd > 0.1 or ystd > 0.1) and (xstd <= prevstds[0] and ystd <= prevstds[1]):
        return plot_image(posx-dx+xavg, posy-dy+yavg, count, [xstd,ystd])
    else:
        return posx-dx+xavg, posy-dy+yavg, 1/(xstd**2), 1/(ystd**2), count

### THINGS FOR RSKY #######

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x - center) / sigma) ** 2 * (-0.5))


# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals) ** 2) / nbin


### THINGS FOR MLR #######
# flux makes magnitude of a given flux in janskys.
def mag(val):
    return -2.5 * np.log10(val / 3631)

##### ACTUAL PROGRAM BITS #####

# define the directory that contains the images
dir = os.environ['HSTDIR']

# setting up arrays with three elements, all zeros - placeholders
#filters = ['F475W', 'F814W', 'F160W']
filters = ['F475W', 'F814W']

data = [0 for x in range(len(filters))]
header = [0 for x in range(len(filters))]
fnu = [0 for x in range(len(filters))]
exp = [0 for x in range(len(filters))]
gain = [0 for x in range(len(filters))]
dark = [0.0022, 0.0022, 0.045]
RN = [0 for x in range(len(filters))]

# width MUST be odd.
width = 15000
pixr = int((width - 1) / 2)
radii = np.arange(40) + 1
#radii = np.arange(30) + 1
area = [0 for x in range(len(radii))]

solarLum = 3.846 * 10 ** 33

# specify the position of the science target and the size of the region around the science target to consider
wavelengths = np.array([475, 814, 1600])  # *u.nm
# wavelengths = np.array([475, 814])

flux = np.zeros([len(filters), len(radii)])  # *u.Jy
subflux = np.zeros([len(filters), len(radii)])
fluxpix = np.zeros([len(filters), width, width])

pixsky = np.zeros([len(filters), width, width])
pixshot = np.zeros([len(filters), width, width])
pixread = np.zeros([len(filters), width, width])
pixdark = np.zeros([len(filters), width, width])
pixNoise = np.zeros([len(filters), width, width])
SNR = np.zeros([len(filters), width, width])

totalphotons = 0

rSky = np.zeros([len(galaxies), len(filters)])

# calculate area of each bagel
for i in range(0, len(area)):
    area[i] = math.pi * math.pow(radii[i], 2)
    if i == 0:
        area[i] = math.pi * math.pow(radii[0], 2)
    else:
        area[i] = math.pi * (math.pow(radii[i], 2) - math.pow(radii[i - 1], 2))

with PdfPages('ad_MONSTER.pdf') as pdf:
    for w in range(0, len(galaxies)):

        # for w in range(0,len(galaxies)):
        # print(galaxies[w].name)
        mxs = [0, 0, 0]
        mys = [0, 0, 0]
        mstdx = [0, 0, 0]
        mstdy = [0, 0, 0]
        for i in range(0, 2):#len(filters)):
            # file = glob.glob(dir+galaxies[w].name+'_final_'+filters[i]+'*sci.fits')
            file = glob.glob(dir + galaxies[w].name + '*/' + res[0] + '/' + str(filters[i]) + '/' + type[1] + '_' + str(
                filters[i]) + '*sci.fits') #res[0]=fine; res[1]=coarse

            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
            RN[i] = header[i]['READNSEA']
            # CENTROID STUFF AKA SG STEAL YO BIIIITCH
            # define positions for photometry
            positions = [galaxies[w].x, galaxies[w].y]

            mxs[i], mys[i], mstdx[i], mstdy[i], count = plot_image(positions[0], positions[1], 0, [100, 100])
        galaxies[w].x = np.average(mxs, weights=mstdx)
        galaxies[w].y = np.average(mys, weights=mstdy)
        # ws.write(w+1, 0, galaxies[w].name)
        # ws.write(w+1, 1, galaxies[w].z)
        ws.write(w + 1, 2, galaxies[w].x)
        ws.write(w + 1, 3, galaxies[w].y)

        print(galaxies[w].x, galaxies[w].y)
        ###    ##I now think i may need to exit and reenter the program? no that's not right.. we could do a different center for each wavelength? i'll ask aleks.
        for i in range(0, len(filters)):
            # GAU IMG BKG TO FIND RSKY TERM:
            # select all pixels with meaningful information and put them into a
            # one-dimensional array
            pixels = data[i][(data[i] == data[i])].ravel()

            # generate a histogram of pixel values for this array
            bin_lo = -100
            bin_hi = 100
            bin_edges = np.linspace(bin_lo, bin_hi, 5000)
            hist, edges = np.histogram(pixels, bins=bin_edges)
            bin_mids = (bin_edges[0:-1] + bin_edges[1:]) / 2.

            # decide on range of pixel values for gaussian fit
            rlo = -5
            rhi = 5
            use = (bin_mids > rlo) & (bin_mids < rhi)
            nbin = len(use)

            # estimate the errors in each bin by assuming poisson statistics: sqrt(n)
            errors = np.sqrt(hist)

            # decide on initial parameters for the gaussian
            params_initial = [1e4, 0, 10]

            # find the best-fit gaussian
            results = minimize(objective_function, params_initial, method='Powell')
            params_fit = results.x
            #model = gaussian(bin_mids, *params_fit)
            rSky[w][i] = params_fit[2]

            fluxpix[i] = data[i]
            pixNoise[i] = np.sqrt((rSky[w][i]) ** 2 + np.abs(fluxpix[i]) + RN[i] ** 2 + dark[i] * exp[i])
            pixshot[i] = np.sqrt(np.abs(fluxpix[i]))
            pixsky[i] = rSky[w][i]
            pixread[i] = RN[i]
            pixdark[i] = np.sqrt(dark[i] * exp[i])
            SNR[i] = fluxpix[i] / pixNoise[i]

            #storing sigma images
            hdu1 = hdu
            hdu1[0].data = pixNoise[i]
            sigma_name = galaxies[w].name + '_' +filters[i] + '_' + res[0] + '_sigma'+'.fits'
            hdu1.writeto(sigma_name)

        # for i in range(0, len(wavelengths)):
        #     # Linear signal to noise ratio
        #     fig = plt.figure()
        #     plt.imshow(SNR[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'SNR_linear' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #
        #     fig = plt.figure()
        #     # plt.imshow(colorUV)
        #     snrimg = img_scale.log(SNR[i], scale_min=1, scale_max=500)
        #     plt.imshow(snrimg)
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'SNR_log' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #     wb.save('updated.xlsx')
        #
        #     # Plotting VJ color map and pixel restrictions.
        #     fig = plt.figure()
        #     plt.imshow(pixNoise[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixNoise_linear' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     noiseimg = img_scale.log(pixNoise[i], scale_min=1, scale_max=500)
        #     plt.imshow(noiseimg)
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixNoise_log' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     # plot the flux values
        #     fig = plt.figure()
        #     plt.imshow(fluxpix[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'fluxpix_linear' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     fluximg = img_scale.log(fluxpix[i], scale_min=1, scale_max=500)
        #     plt.imshow(fluximg)
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'fluxpix_log' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     plt.imshow(pixshot[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixshot' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     shotimg = img_scale.log(pixshot[i], scale_min=1, scale_max=500)
        #     plt.imshow(shotimg)
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixshot' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     plt.imshow(pixsky[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixsky' + str(i))
        #     pdf.savefig()
        #     plt.close()
        #
        #     fig = plt.figure()
        #     plt.imshow(pixdark[i])
        #     plt.colorbar()
        #     plt.title(galaxies[w].name + 'pixdark' + str(i))
        #     pdf.savefig()
        #     plt.close()

os.system('open %s &' % 'ad_MONSTER.pdf')