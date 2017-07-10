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

# flux makes magnitude of a given flux in janskys.
def mag(val):
    return -2.5*np.log10(val/3631)

# flux makes luminosity defined by Hogg 1999??
def luminosity(wavelength,flux,lDcm):
    return const.c*u.s/u.m/(wavelength*10**-9)*flux*10**-23*(4*math.pi*lDcm**2)/solarLum

# Given two coefficients and the color will return the mass-to-light ratio
# as defined by Bell and DeJong 2000
def mLR(a,b,color):
    return 10**(a+b*color)


# define the directory that contains the images
dir = os.environ['HSTDIR']
#dir = '/Users/aaw/data/'

Ba = [-1.019,-1.113,-1.026,-.990,-1.110,-.994,-.888]
Bb = [1.937,2.065,1.954,1.883,2.018,1.804,1.758]
B_coeff = [Ba,Bb]
Va = [-.759,-.853,-.766,-.730,-.850,-.734,-.628]
Vb = [1.537,1.665,1.554,1.483,1.618,1.404,1.358]
V_coeff = [Va,Vb]
Ja = [-.540,-.658,-.527,-.514,-.659,-.621,-.550]
Jb = [.757,.907,.741,.704,.878,.794,.801]
J_coeff = [Ja,Jb]
coeff = [B_coeff, V_coeff, J_coeff]


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

solarLum = 3.846*10**33

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
lDcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc
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
    #for w in range(0,1):
    for w in range(0,len(galaxies)):
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
        n = np.ma.masked_where(SNR<10,fluxpix)
        for i in range(0,len(filters)):
            n[i]=n[i]*fnu[i]/exp[i]

        colorUV = mag(n[0])-mag(n[1])
        for i in range (0,len(filters)):
            lum = luminosity(filters[i], n, lDcm[w])
            for mod in range(5,6):
                fig = plt.figure()
                mlr = mLR(coeff[i][0][mod], coeff[i][1][mod], colorUV)
                mass = lum*mlr
                plt.imshow(m[i]*mass[i])
                plt.colorbar()
                plt.title(galaxies[w]+' Mass Profile in Solar Masses at ' + wavelengths[i], fontsize = 12)
                pdf.savefig()
                plt.close()
                            
        fig = plt.figure()
        plt.imshow(colorUV)
        plt.colorbar()
        plt.title('ColorUV map and constrictions')
        pdf.savefig()
        plt.close()
