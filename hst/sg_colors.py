""" 
Sophia C W Gottlieb I
sg_colors.py
20170808

FOR EACH GALAXY colors will plot comparisons of U-V and V-J colors
using import statements of sg_monster, SNR, 
We will NOT be using centroid in this particular code.
"""
import os
import numpy as np
from astropy.io import fits
import glob
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
from scipy import integrate
from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

# define the directory that contains the images
dir = os.environ['HSTDIR']
#dir = '/Users/sgottlie/Desktop/'
#dir = '/Volumes/physics/linux-lab/data/hst/'

### THINGS FOR READING IN GALAXY DATA FROM galaxydata.xlsx

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

from mmap import mmap,ACCESS_READ
from xlrd import open_workbook

# for writing into the xl file
import xlwt
"""
#Define Galaxy class to hold name, redshift, x and y positions for fine and
coarse data, and have functions for luminosity distance in cms and rad to 
kpc conversion 
Okay so i checked and the functions don't work and I hate them. 
Anyway, I saved them as fields and it's working fine now
"""
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
cf = ['coarse','fine']

# label describing image processing
rc = ['final','convolved_image']
# matching file name ending whoooooooot.
# not confident the h in needed in whooooot. but lets roll with it for now. its a nice h
fex = ['*sci.fits','.fits']

#wbo = open_workbook(res[t]+'data.xlsx')
wbo = open_workbook('coarsedata.xlsx')
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

####THINGS FOR COLOR
# flux makes magnitude of a given flux in janskys.
def mag(val):
    return -2.5*np.log10(val/3631)

titleMFKRs = ['UV coarse raw','UV coarse conv','UV fine raw','UV fine conv','VJ coarse raw','VJ coarse conv']
## THINGS FOR SNR
# RSKY:

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

#setting up arrays with three elements, all zeros - placeholders
filters = ['F475W','F814W','F160W']
dark = [0.0022,0.0022,0.045]

# Setting up the bounds of the analysis area:
# width MUST be odd.
width = 81
pixr = int((width-1)/2)
#we have decided to stack things??? does it actually matter??? recent studies suggest no. lets push a work around for our ridikidonk friend 1600nm and cruise forward with the linear data structure keeping in mind that index of 7 will be bare and 1600 c f will be in index 8 aka the structure will be 9 in length
fluxpix = np.zeros([9, width, width])
# I am not 100% convinced I am going to need 9 of these... I think we only need that for fluxpix...

pixNoise = np.zeros([width,width])
SNR = np.zeros([width,width])

#Here ye, here ye. MY name is Samuel Sabrey and I present pairs of data files to be compared to form colors we are interested in:
pairs = [[0,4],[1,5],[2,6],[3,6],[5,8],[6,8]]

#THIS MIGHT NOT WORK IDK THO. ALSO UNSURE OF HOW THE PAIRS WITHIN PARENTHESES WORK LOL
def color (pair):
    return mag(fluxpix[pair[0]])-mag(fluxpix[pair[1]])

with PdfPages('sg_COLORS.pdf') as pdf:
    #for w in range(0,1):
    for w in range(0, len(galaxies)):
        J = galaxies[w]
        print(J.name)
        #for i in range(0,1):
        for i in range (0,len(filters)):
            for res in range(0,2):
                for prc in range (0,2):
                    ih8u = i*4+res*2+prc
                    print(ih8u)
                    # i would love to through a 'does this file exist' thinggy in here.
                    #nameYoFile = dir+J.name+'/'+cf[res]+'/'+filters[i]+'/'+rc[prc]+'_'+filters[i]+fex[prc]
                    #file = glob.glob(nameYoFile)
                    if ih8u != 7 and ih8u != 9 and ih8u != 10 and ih8u != 11 and ih8u != 12:
                        file = glob.glob(dir+J.name+'/'+cf[res]+'/'+filters[i]+'/'+rc[prc]+'_'+filters[i]+fex[prc])

                        hdu = fits.open(file[0])
                        data, header = hdu[0].data, hdu[0].header
                        fnu = header['PHOTFNU']
                        exp = header['EXPTIME']
                        gain = header['CCDGAIN']
                        RN = header['READNSEA']

                        """RSKY THINGS: happens for each data set, no need to hold onto used data.

    select all pixels with meaningful information and put them into a one-dimensional array """
                        pixels = data[(data==data)].ravel()

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
                        rSky = params_fit[2]
                        ###

                        ### SNR code: we loop through every pixel because of how things get tilted in the code.
                        for j in range(0,width):
                            for k in range(0,width):
                                # Here's a little story. It's hard to figure out what center we need because I set it up dumbly. LUCKILY the fine center is LITERALLY the coarse center *2. The way I manage that is by multiplying the coarse center by res+1. for res+1, coarse(0) is 1 and fine(1) is two. VOILAAAA!!! now I have to go punch my code in the face upwards
                                #ALSO please remember that these indices are SUPER FUCKED UP because of the difference in indexing between DS9 and python. I think we could take or leave the -1s, that may have just been me, fucking some shit up.... my bad... we should talk about that because it feels wrong.....
                                #this feels like something I could have used a stamp for. hush... later....

                                fluxpix[ih8u,width-j-1,k] = data[int(j+round(J.y*(res+1))-pixr+1)][int(k+round(J.x*(res+1))-pixr+1)]
                                #pixnoise may not need to be the same structure as fluxpix. infact, let us attempt that:
                                # OR: pixNoise[ih8u,width-j-1,k] = math.sqrt((rSky)**2+fluxpix[ih8u][width-j-1][k]+(RN**2+(gain/2)**2)*1+dark[i]*exp)
                                val = (rSky)**2+fluxpix[ih8u][width-j-1][k]+(RN**2+(gain/2)**2)*1+dark[i]*exp
                                if val > 0:
                                    pixNoise[width-j-1,k] = math.sqrt((rSky)**2+fluxpix[ih8u][width-j-1][k]+(RN**2+(gain/2)**2)*1+dark[i]*exp) #the matching indices keeps things in the correct shapes
                                # Incorrect because it doesn't map properly: pixNoise[j,k] = math.sqrt((rSky)**2+fluxpix[ih8u][j][k]+(RN**2+(gain/2)**2)*1+dark[i]*exp)
                                #SNR[ih8u, width-j-1,k] = fluxpix[ih8u,width-j-1,k]/pixNoise[ih8u,width-j-1,k]
                                    SNR[width-j-1,k] = fluxpix[ih8u,width-j-1,k] / pixNoise[width-j-1,k]
                                else:
                                    SNR[width-j-1,k] = -10    

                                fluxpix[ih8u] = np.ma.masked_where(SNR<3, fluxpix[ih8u])

                            #i think that completes the data collection, snr, mask, store & exit stage.... that is done for each galaxy. now we return to the galaxy level.
        # and i believe that is here. so at this point, we have our fluxpix that is 0-6 and 8. meow.
        fig = plt.figure()
        plt.suptitle(J.name +' fliggity flux, dawg')
        for h0t in range(0,len(fluxpix)):
            ax = fig.add_subplot(3,3,h0t+1)
            plt.imshow(fluxpix[h0t])
            #plt.tight_layout()
            plt.colorbar(format='%.0e')
            f = int(h0t/4)
            r = int((h0t-f*4)/2)
            im = int(h0t%2)
            
            plt.title(filters[f]+' '+cf[r]+ ' ' +rc[im])
        pdf.savefig()
        plt.close() #??
            #plt.close()
            #idt i want to close the plt just yet. i think that makes a whole new thingy and i'm too young for that shit tbh. uncertain if pdf.savefig() is needed either. let's take it down a level for now
        fig = plt.figure()
        plt.suptitle(J.name +' Color Comparison')
        for p00ps in range(0,len(pairs)):
            ax = fig.add_subplot(2,3,p00ps+1)
            pl0tme = color(pairs[p00ps])
            plt.imshow(pl0tme) #remember how to make red / blue!!!!
            plt.colorbar()
            plt.tight_layout()
            ax.set_title(titleMFKRs[p00ps]) #SRY!SRY
        pdf.savefig()
        plt.close() #???
#I wouldn't hate it if I at some point learned to make titles that weren't like, WILDLY inappropriate.
            
        
