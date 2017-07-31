# Sophia Gottlieb
# 20170710

# sg_monster.py
# This code will pull out all the stops
# If possible, it will begin with centroid code.
# Then, SNR for pixel
# Then, MLR for significant pixels
# But to begin, we will do SNR and MLR because those are mine....
# currently does SNR and MLR in pixel mode....
# need to talk to aleks about SNR restrictions via filter vs color....
# Also now includes GAU code for rSky term.
# Up next we will include the cl_centroidmerged.py
# It now reads in a file, and stores them as Galaxies
# Next i have to make it do annular stuff.
# I guess I am now trying to incorporate the centroid stuff into 
#import relavent packages
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
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc
# setting up final xl file
wb = xlwt.Workbook()
ws = wb.add_sheet('Galaxy Data')

# Grabbing and filling galaxy data
res = ['fine','coarse']
t=0

#type = ['final','final','final']
#type = ['final','final']
#type = ['convolved_image', 'convolved_image', 'final']
type = ['convolved_image', 'final']

#wbo = open_workbook(res[t]+'data.xlsx')
wbo = open_workbook('galaxydata.xlsx')
#wbo = open_workbook('update.xls')
for sheet in wbo.sheets():
    numr = sheet.nrows
    galaxies = []
    
    ws.write(0, 0, 'name')
    ws.write(0, 1, 'z')
    ws.write(0, 2, 'x')
    ws.write(0, 3, 'y')
    
    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))
        ws.write(row,0,sheet.cell(row,0).value)
        ws.write(row,1,sheet.cell(row,1).value)


###
### THINGS FOR CENTROID ###
dx = 5
dy = dx

# Maybe we could make a weighted average? Idk.
# plot_image is our centroid code. I have defined it to be recursive, but it
# seems that we don't need it. we give plot image the estimated position in x and y
# in relation to the whole image
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
    print(count, posx-dx+xavg, posy-dy+yavg, xstd, ystd)
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
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin


### THINGS FOR MLR #######
# flux makes magnitude of a given flux in janskys.
def mag(val):
    return -2.5*np.log10(val/3631)

# flux makes luminosity defined by Hogg 1999??
def luminosity(wavelength,flux,lDistcm):
    return const.c*u.s/u.m/(wavelength*10**-9)*flux*10**-23*(4*math.pi*lDistcm**2)/solarLum

# Given two coefficients and the color will return the mass-to-light ratio
# as defined by Bell and DeJong 2000
def mLR(a,b,color):
    return 10**(a+b*color)
##### ACTUAL PROGRAM BITS #####

# define the directory that contains the images
dir = os.environ['HSTDIR']


# coefficients from Bell & de Jong 2001, note: our photometry does not correspond exactly with restframe B, V, and J.
# We could apply k-corrections to estimate restframe magnitudes at these wavelengths but are not doing so currently. 
# The plan going forward is to use uhm... models that do not require k-corrections. 
# Don't uhm me. Aleks is the best. Those models are iSEDfit ft. Moustakas and Prospector. 
Ba = [-1.019,-1.113,-1.026,-.990,-1.110,-.994,-.888]
Bb = [1.937,2.065,1.954,1.883,2.018,1.804,1.758]
B_coeff = [Ba,Bb]
Va = [-.759,-.853,-.766,-.730,-.850,-.734,-.628]
Vb = [1.537,1.665,1.554,1.483,1.618,1.404,1.358]
V_coeff = [Va,Vb]
Ja = [-.540,-.658,-.527,-.514,-.659,-.621,-.550]
Jb = [.757,.907,.741,.704,.878,.794,.801]
J_coeff = [Ja,Jb]
# coeff has three bracket sets: [wavelength][a coef][b coef]
coeff = [B_coeff, V_coeff, J_coeff]

#setting up arrays with three elements, all zeros - placeholders
filters = ['F475W','F814W','F160W']
#filters = ['F475W','F814W']

data = [0 for x in range(len(filters))]
header = [0 for x in range(len(filters))]
fnu = [0 for x in range(len(filters))]
exp = [0 for x in range(len(filters))]
gain = [0 for x in range(len(filters))]
dark = [0.0022,0.0022,0.045]
RN = [0 for x in range(len(filters))]

#width MUST be odd.
width = 81
pixr = int((width-1)/2)
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]
area = [0 for x in range(len(radii))]
annNoise = np.zeros([len(filters),len(radii)])
annSNR = np.zeros([len(filters),len(radii)])


solarLum = 3.846*10**33

# specify the position of the science target and the size of the region around the science target to consider
wavelengths = np.array([475, 814, 1600]) #*u.nm
#wavelengths = np.array([475, 814])

flux = np.zeros([len(filters),len(radii)]) #*u.Jy
subflux = np.zeros([len(filters),len(radii)])
fluxpix = np.zeros([len(filters), width, width])

pixNoise = np.zeros([len(filters), width, width])
SNR = np.zeros([len(filters), width, width])

totalphotons = 0

rSky = np.zeros([len(galaxies),len(filters)])

#calculate area of each bagel
for i in range(0, len(area)):
    area[i] = math.pi*math.pow(radii[i],2)
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))


with PdfPages('sg_MONSTER.pdf') as pdf:
    #for w in range(0,1):

    for w in range(0,len(galaxies)):
        print(galaxies[w].name)
        mxs = [0,0,0]
        mys = [0,0,0]
        mstdx = [0,0,0]
        mstdy = [0,0,0]
        for i in range(0, len(filters)):
            #file = glob.glob(dir+galaxies[w].name+'_final_'+filters[i]+'*sci.fits')
            file = glob.glob(dir+galaxies[w].name+'*/'+res[1]+'/'+str(filters[i])+'/'+type[1]+'_'+str(filters[i])+'*sci.fits')

            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
            RN[i] = header[i]['READNSEA']
# CENTROID STUFF AKA SG STEAL YO BIIIITCH
            #define positions for photometry
            positions = [galaxies[w].x, galaxies[w].y]
            
            mxs[i], mys[i], mstdx[i], mstdy[i], count = plot_image(positions[0], positions[1], 0, [100,100])
        galaxies[w].x = np.average(mxs, weights = mstdx)
        galaxies[w].y = np.average(mys, weights = mstdy)
        #ws.write(w+1, 0, galaxies[w].name)
        #ws.write(w+1, 1, galaxies[w].z)
        ws.write(w+1, 2, galaxies[w].x)
        ws.write(w+1, 3, galaxies[w].y)

        print(galaxies[w].x, galaxies[w].y)
        ###    ##I now think i may need to exit and reenter the program? no that's not right.. we could do a different center for each wavelength? i'll ask aleks.
        for i in range(0,len(filters)):
# GAU IMG BKG TO FIND RSKY TERM: 
            # select all pixels with meaningful information and put them into a
            # one-dimensional array
            pixels = data[i][(data[i]==data[i])].ravel()

            # generate a histogram of pixel values for this array
            bin_lo = -100
            bin_hi = 100
            bin_edges = np.linspace(bin_lo, bin_hi, 5000)
            hist, edges = np.histogram(pixels, bins=bin_edges)
            bin_mids = (bin_edges[0:-1] +  bin_edges[1:]) / 2.
    
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
            model = gaussian(bin_mids, *params_fit)
            rSky[w][i]= params_fit[2]
# BEGINNING SNR CODE: 
            
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    # In data numbers? In electrons?
                    fluxpix[i,width-j-1,k] = data[i][int(j+round(galaxies[w].y)-pixr+1)][int(k+round(galaxies[w].x)-pixr+1)]
                    pixNoise[i,width-j-1,k] =  math.sqrt((rSky[w][i]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i])
                                              
                    SNR[i,width-j-1,k] = fluxpix[i,width-j-1,k]/pixNoise[i,width-j-1,k]

            #do photometry on images
            #convert to proper units
            for j in range(0,len(radii)):
                aperture = CircularAperture([galaxies[w].x, galaxies[w].y], radii[j])
                phot_table = aperture_photometry(data[i], aperture)
                #flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])/(3.631*10**-6)
                flux[i,j] = phot_table['aperture_sum'][0]#*gain[i]*exp[i]
                if j == 0:
                    subflux[i,j] = flux[i,j]
                else:
                    subflux[i,j] = flux[i,j]-flux[i,j-1]
        
                annNoise[i,j] = np.sqrt((rSky[w][i]*area[j])**2+flux[i][j]+(RN[i]**2)*area[j]+dark[i]*area[j]*exp[i])
                annSNR[i,j] = flux[i,j]/annNoise[i,j]

        m = np.ma.masked_where(SNR<10,SNR)
        n = np.ma.masked_where(SNR<10,fluxpix)
# BEGINNING MLR CODE BUT ONLY IN PIXELS???: 
        for i in range(0,len(wavelengths)):
            m[i] = np.ma.masked_where(SNR[0]<10,SNR[i])
            n[i] = np.ma.masked_where(SNR[0]<10, fluxpix[i])
            n[i]=n[i]*fnu[i]/exp[i]
        # making uv color
        colorUV = mag(n[0])-mag(n[1])
        for i in range (0,len(wavelengths)):
            lD = galaxies[w].lDcm
            lum = luminosity(wavelengths[i], n, lD)
            # model is fixed to star formation epoch with bursts
            for mod in range(5,6):
                fig = plt.figure()
                mlr = mLR(coeff[i][0][mod], coeff[i][1][mod], colorUV)
                mass = lum*mlr
                #masked colorplot of mass. 
                plt.imshow(m[i]*mass[i])
                #plt.imshow(mass)
                plt.colorbar()
                plt.title(galaxies[w].name+' Mass Profile in Solar Masses at ' + filters[i], fontsize = 12)
                pdf.savefig()
                plt.close()
        # Plotting UV color map and pixel restrictions.                     
        fig = plt.figure()
        plt.imshow(colorUV)
        plt.colorbar()
        plt.title(galaxies[w].name+'ColorUV map and constrictions')
        pdf.savefig()
        plt.close()
        wb.save('updated.xlsx')
        
        colorVJ = mag(n[1])-mag(n[2])
        # Plotting VJ color map and pixel restrictions.                     
        fig = plt.figure()
        plt.imshow(colorVJ)
        plt.colorbar()
        plt.title(galaxies[w].name+'ColorVJ map and constrictions')
        pdf.savefig()
        plt.close()
        
os.system('open %s &' % 'sg_MONSTER.pdf')
