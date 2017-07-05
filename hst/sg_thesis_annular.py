# SCGWI 20170623
# Just trying to cut down on unnecessary code
#

# import modules of Python
import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from scipy import integrate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

# set up constants for the program
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
xcen = [3628, 3933, 3386.5, 3477.5, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycen = [4153, 4136, 3503.2, 3404.3, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]
zs = [0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752]
lDcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc

#collection = ['F475W','F814W']
#filters = np.array([475, 814])

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

radToKpc = conv.arcsec_per_kpc_proper(zs)*0.05/u.arcsec*u.kpc

collection = ['F475W','F814W','F160W']
filters = np.array([475, 814, 1600])

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

# specifics for the program
solarLum = 3.846*10**33    #solar Mass is the mass of the sun in kg
#radToKpc = 0.05*7.194      #converts radius to kpc

radii = np.arange(40)+1
#radii = radii * radToKpc
area = [0 for x in range(len(radii))]

flux = np.zeros([len(collection),len(radii)]) #*u.Jy
subflux = np.zeros([len(collection),len(radii)])

colorUV, colorVJ = np.zeros([len(collection)-1, len(radii)])
colors = [colorUV, colorVJ]
#colorUV = np.zeros([len(collection)-1, len(radii)])
#colors = [colorUV]

#flux = np.zeros([len(radii)])
#subflux = np.zeros([len(radii)])
dir = os.environ['HSTDIR']

# set up lots of data structures
data = [0 for x in range(len(filters))]
header = [0 for x in range(len(filters))]
fnu = [0 for x in range(len(filters))]
exp = [0 for x in range(len(filters))]
#gain = [0 for x in range(len(filters))]

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

pltcolors = ['b','g','r','c','m','y','k']
pltmarker = ['o','*','v','P','h','d']
# loop filters
# grab data
# We pass in file name, positions, and does photometry and returns subflux
def readDataGrabFlux(fileName, positions):
    hdu = fits.open(file[0])
    data[i], header[i] = hdu[0].data, hdu[0].header
    fnu[i] = header[i]['PHOTFNU']
    exp[i] = header[i]['EXPTIME']
    for j in range(0,len(radii)):
        aperture = CircularAperture(positions, radii[j])
        phot_table = aperture_photometry(data[i], aperture)
        flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])
        if j == 0:
            subflux[i,j] = flux[i,j]
        else:
            subflux[i,j] = flux[i,j]-flux[i,j-1]
    return flux, subflux

def totalFlux(flux):
    return flux[len(radii)-1]
# flux makes mag and color -> MLR

def mag(val):
    return -2.5*np.log10(val/3631)
# flux makes luminosity

def luminosity(wavelength,flux,lDcm):
    return const.c*u.s/u.m/(wavelength*10**-9)*flux*10**-23*(4*math.pi*lDcm**2)/solarLum

def mLR(a,b,color):
    return 10**(a+b*color)

# luminosity MLR makes mass.
def mass(mLR, luminosity):
    return mLR*luminosity
#mass(mLR(a,b,color), luminosity(wavelength[i],flux, lDcm)
# graph some stuff
    

# decide to plot things
with PdfPages('sg_THESIS.pdf') as pdf: 
    # calculation of mass by galaxy
    for gal in range(0,len(galaxies)):
        kpcArea = area*radToKpc[gal]**2
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.title('1000s of Mass Density vs Radius for ' + galaxies[gal])
        plt.xlabel('Radii (kpc)')
        plt.ylabel('Mass (solar mass / kpc^2)')
        # grab all flux data by looping through wavelengths
        for i in range(0,len(collection)):
            file = glob.glob(dir+galaxies[gal]+'_final_'+collection[i]+'*sci.fits')
            flux, subflux = readDataGrabFlux(file,[xcen[gal],ycen[gal]])
        # grab all radii to get magnitude / color
        for r in range(0,len(radii)):
            # make both colors
            for c in range(0,len(colors)):
                colors[c][r] = mag(subflux[c,r])-mag(subflux[c+1,r])
        # analysis for all filters
        for i in range(0,len(filters)):
        #for i in range(1,2):
            lum = luminosity(filters[i], subflux,lDcm[gal])
            # for both colors
            #for c in range(0,1):
            for c in range(0,len(colors)):
                # for all seven models of B&dJ
                #for m in range(0,7):
                for m in range(5,6):
                    desc = 'f' + str(i) + 'c' + str(c) + 'm' + str(m)
                    mlr = mLR(coeff[i][0][m],coeff[i][1][m], colors[c])
                    mass = lum*mlr
                    ax.plot(radii*radToKpc[gal], mass[i]/kpcArea,linestyle = 'None',  marker = pltmarker[c], color = pltcolors[i], label = desc)
                    legent = ax.legend(loc='upper right')
        pdf.savefig()
        plt.close()
                #print('mass calculation')
        # Up next we need to figure out how to graph all 102010483857205203948 calcualtions we made of the mass. 

