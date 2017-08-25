""" Sophia C W Gottlieb I
my attempt to understand binning. so we are just going to grab a single
data file so it doesn't take me a million years to run edits...
Also TBH i'm deciding whether or not this is going to steal a huge
chunk of SNR code because IRL that's what I will be working with...
Might as well, right? k. fuck.
"""


#20170825

import os
dir = os.environ['HSTDIR']

import numpy as np
#fuck i really don't want to copy paste all this shit. fuck this.

from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as cost
from astropy.cosmology import FlatLambdaCDM
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
#We define a galaxy data class so I don't have to fun
class Galaxy:
    def __init__(self, name, z, x,y):
        self.z = z
        self.y = y
        self.x = x
        self.name = name
        #lDcm gives the luminosity distance in centimeters
        #self.lDcm = cosmo.luminosity_distance(0.603)*u.Mpc.to(u.cm)/u.Mpc
        #self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        #pixtoKpc is the conversion betweem pixels to rad to kpc based on z value
        #self.pixtoKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc
        self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc
# Let's set up some stuff so that I can grab stuff from inside an excel sheet.

from xlrd import open_workbook
wbo = open_workbook('galaxydata.xlsx')
for sheet in wbo.sheets():
    numr = sheet.nrows
    galaxies = []

    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))

# i don't think we need centroid stuff still... we still need rsky stuff tho :P

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

## have some stuff!
filters = ['F475W','F814W','F160W']
wavelengths = np.array([475, 814, 1600])
dark = [0.0022,0.0022,0.045]
data = [0 for x in range(len(filters))]
header = [0 for x in range(len(filters))]
fnu = [0 for x in range(len(filters))]
exp = [0 for x in range(len(filters))]
gain = [0 for x in range(len(filters))]

RN = [0 for x in range(len(filters))]

#width MUST be odd.
width = 81
pixr = int((width-1)/2)
radii = np.arange(40)+1

#the fuck is this????? do i need it???
solarLum = 3.846*10**33

fluxpix = np.zeros([len(filters), width, width])

pixNoise = np.zeros([len(filters), width, width])
SNR = np.zeros([len(filters), width, width])

totalphotons = 0

#rSky = np.zeros([len(galaxies),len(filters)])
rSky = [[ 10.31036745,  10.64816775,   7.50788545],
       [ 10.5783909 ,  13.35756601,   7.58130161],
       [ 11.77636022,  10.54219567,   8.27527023],
       [ 11.28425007,  12.25892311,   9.03878993],
       [ 11.90782971,  13.34688462,   8.67763722],
       [ 10.62194407,  10.75379768,   7.62254201],
       [ 10.62821474,  10.88078454,   7.70920673],
       [ 10.94261798,  10.26098238,   6.74143006],
       [ 11.35125336,   9.69418394,   6.87375846],
       [ 10.5804402 ,  10.592764  ,   7.46983987],
       [  9.08922928,  10.57481831,   7.83102976],
       [ 11.18114582,   8.79336163,   8.10811507]]

import glob
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
res = 0

with PdfPages('sg_BINS.pdf') as pdf:
    for w in range(0,1):
        J = galaxies[w]
        for i in range(0,1):
        #for i in range(0,len(filters)):
            file = glob.glob(dir+J.name+'/coarse/'+filters[i]+'/final_'+filters[i]+'*sci.fits')
            
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
            fnu = header['PHOTFNU']
            exp = header['EXPTIME']
            gain = header['CCDGAIN']
            RN = header['READNSEA']

            #onto gau_img to find rsky. :P should i just import this? for now? idk. yeah
            #rSky = 10.31036745
            # SNR CODE!
            for j in range(0,width):
                for k in range(0,width):

                    fluxpix[i,width-j-1,k] = data[int(j+round(J.y*(res+1))-pixr+1)][int(k+round(J.x*(res+1))-pixr+1)]
                    val = ([w][i])**2+fluxpix[i][width-j-1][k]+(RN**2+(gain/2)**2)*1+dark[i]*exp
                    if val > 0:
                        pixNoise[i,width-j-1,k] = np.sqrt(val) 
                        SNR[i,width-j-1,k] = fluxpix[i,width-j-1,k] / pixNoise[i,width-j-1,k]
                    else:
                        SNR[i,width-j-1,k] = -10    

                    fluxpix[i] = np.ma.masked_where(SNR[i]<3, fluxpix[i])
