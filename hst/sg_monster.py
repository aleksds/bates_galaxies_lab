# Sophia Gottlieb
# 20170710

# sg_monster.py
# This code will pull out all the stops
# If possible, it will begin with centroid code.
# Then, SNR for pixel
# Then, MLR for significant pixels
# But to begin, we will do SNR and MLR because those are mine....

#import relavent packages
import os
import numpy as np
from astropy.io import fits
#from photutils import CircularAperture
#from photutils import aperture_photometry
import glob
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
#dir = os.environ['HSTDIR']
dir = '/Users/aaw/data/'

#setting up arrays with three elements, all zeros - placeholders
wavelengths = ['F475W','F814W','F160W']

data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]
gain = [0 for x in range(len(wavelengths))]
dark = [0.0022,0.0022,0.045]
RN = [0 for x in range(len(wavelengths))]

#width MUST be odd.
width = 81
pixr = int((width-1)/2)
rSky = np.zeros([len(wavelengths),12])
rSky[0]= (10.31036745,  10.5783909 ,  11.77636022,  11.28425007,11.90782971,  10.62194407,  10.62821474,  10.94261798, 11.35125336,  10.58044017, 9.08922928,  11.18114582)
rSky[1]= (10.64816775,  13.35756601,  10.54219567,  12.25892311,13.34688468,  10.75379768,  10.88078454,  10.26098238, 9.69418394,  10.592764, 10.57481831, 8.79336163)
rSky[2]= (7.50788545,  7.58130161,  8.27527023,  9.03878993,  8.67763722, 7.62254201,  7.70920672,  6.74143006, 6.87375846,  7.46983987, 7.83102976,  8.10811507)


# specify the position of the science target and the size of the region around the science target to consider
filters = np.array([475, 814, 1600]) #*u.nm
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
xcen = [3628, 3933, int(3386.5), 3477, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycen = [4153, 4136, int(3503.2), 3404, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]
zs = [0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752]

#flux = np.zeros([len(wavelengths),len(radii)]) #*u.Jy
#subflux = np.zeros([len(wavelengths),len(radii)])
fluxpix = np.zeros([len(wavelengths), width, width])

pixNoise = np.zeros([len(wavelengths), width, width])
SNR = np.zeros([len(wavelengths), width, width])

#Ldcm is the luminosity distance in cm, even though astropy thinks it is in Mpc. 
Ldcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
radToKpc = conv.arcsec_per_kpc_proper(zs)*0.05/u.arcsec*u.kpc
totalphotons = 0


#calculate area of each bagel
#for i in range(0, len(area)):
 #   area[i] = math.pi*math.pow(radii[i],2)
    #if i == 0:
    #    area[i] = math.pi*math.pow(radii[0],2)
    #else:
    #    area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

with PdfPages('sg_MONSTER.pdf') as pdf:
    for w in range(0,1):
        # for w in range(0,len(galaxies)):
        for i in range(0, len(wavelengths)):
            file = glob.glob(dir+galaxies[w]+'_final_'+wavelengths[i]+'*sci.fits')
            #print(file)
            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
            RN[i] = header[i]['READNSEA']

            #define positions for photometry
            positions = [(xcen[w], ycen[w])]
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    # In data numbers? In electrons?
                    fluxpix[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))#*gain[i]*exp[i])
                    SNR[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[i][w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                    pixNoise[i,width-j-1,k] =  math.sqrt((rSky[i][w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i])


        m = np.ma.masked_where(SNR<10,SNR)
        for i in range(0,len(filters)):
            fig = plt.figure()
            plt.imshow(m[i])
            plt.colorbar()
            plt.title(galaxies[w]+' SNR at '+ str(wavelengths[i]))
            pdf.savefig()
            plt.close()
