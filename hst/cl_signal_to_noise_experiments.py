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
from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
import matplotlib.gridspec as gridspec

# define the directory that contains the images
dir = os.environ['HSTDIR']

#setting up arrays with three elements, all zeros - placeholders
wavelengths = ['F475W','F814W','F160W']

data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]
gain = [0 for x in range(len(wavelengths))]
dark = [0.0022,0.0022,0.045]
RN = [0 for x in range(len(wavelengths))]

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
xcen = [3629, 3934, 3387, 3477, 3573, 3803, 3884, 4147, 3787, 4175, 3567, 4067]
ycen = [4154, 4137, 3503, 3405, 3339, 4170, 4165, 3922, 4186, 3827, 3436, 4054]

#from cl_centroidmerged.py
bestxswhole = np.array([ 3628.8630965 ,  3933.65237012,  3386.54069503,  3477.20367975, 3572.97219824,  3803.17972872,  3883.20568985,  4147.00201888, 3786.77909178,  4175.294108  ,  3565.96733481,  4065.88540274])
bestyswhole = np.array([ 4153.8967805 ,  4136.5722187 ,  3503.1663102 ,  3405.39675221, 3339.1874928 ,  4170.32537204,  4164.77089952,  3922.116162  , 4186.05759774,  3826.62884958,  3434.92868183,  4053.58473705])

xcendiff = bestxswhole-xcen
ycendiff = bestyswhole-ycen

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
xcen = [3629, 3934, 3387, 3477, 3573, 3803, 3884, 4147, 3787, 4175, 3567, 4067]
ycen = [4154, 4137, 3503, 3405, 3339, 4170, 4165, 3922, 4186, 3827, 3436, 4054]
zs = [0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752]

# define the radii to be used for aperture photometry
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]
sigskyy = np.zeros([len(wavelengths),len(radii)])
sigread = np.zeros([len(wavelengths),len(radii)])
sigdark = np.zeros([len(wavelengths),len(radii)])
sigstar = np.zeros([len(wavelengths),len(radii)])

flux = np.zeros([len(wavelengths),len(radii)]) #*u.Jy
subflux = np.zeros([len(wavelengths),len(radii)])
fluxpix = np.zeros([len(galaxies),len(wavelengths), width, width])

#fluxvalues = [[0 for x in range(len(wavelengths))] for y in range(len(galaxies))]
pixNoise = np.zeros([len(galaxies),len(wavelengths), width, width])
SNR = np.zeros([len(galaxies),len(wavelengths), width, width])
annNoise = np.zeros([len(wavelengths),len(radii)])

#Ldcm is the luminosity distance in cm, even though astropy thinks it is in Mpc. 
Ldcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
pixtokpc = conv.arcsec_per_kpc_proper(zs)*0.05/u.arcsec*u.kpc
totalphotons = 0

# percunc specifies the percent of variation we expect from systematic error... 
# For now, we have chosen 0.05, or 5%.
percUnc = 0.05

#making list of radius of each pixel
disfrombest=np.zeros([len(galaxies),(width),(width)])
for w in range(0,12):
    for q in range(1,width+1):
        for p in range(1,width+1):
            if p < 41 and q < 41:
                disfrombest[w][p-1][q-1] = ((41-p+xcendiff[w])**2 + (41-q-ycendiff[w])**2)**0.5
            if p == 41 and q < 41:
                disfrombest[w][p-1][q-1] = ((41-p+np.abs(xcendiff[w]))**2 + (41-q-ycendiff[w])**2)**0.5
            if p > 41 and q < 41:
                disfrombest[w][p-1][q-1] = ((41-p-xcendiff[w])**2 + (41-q-ycendiff[w])**2)**0.5
            if p < 41 and q == 41:
                disfrombest[w][p-1][q-1] = ((41-p+xcendiff[w])**2 + (41-q+np.abs(ycendiff[w]))**2)**0.5
            if p == 41 and q == 41:
                disfrombest[w][p-1][q-1] = ((41-p+np.abs(xcendiff[w]))**2 + (41-q+np.abs(ycendiff[w]))**2)**0.5
            if p > 41 and q == 41:
                disfrombest[w][p-1][q-1] = ((41-p-xcendiff[w])**2 + (41-q+np.abs(ycendiff[w]))**2)**0.5
            if p < 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p+xcendiff[w])**2 + (41-q+ycendiff[w])**2)**0.5
            if p == 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p+np.abs(xcendiff[w]))**2 + (41-q+ycendiff[w])**2)**0.5
            if p > 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p-xcendiff[w])**2 + (41-q+ycendiff[w])**2)**0.5
                
# Now, we loop through all galaxies
# adapted for J0905 475 coarse regular flux and residual flux images

for w in range (0, len(galaxies)):
    for i in range (0, 3):
            
        # read in the images
        file = glob.glob(dir+galaxies[w]+'*/coarse/'+wavelengths[i]+'/final*sci.fits')
        hdu = fits.open(file[0]) #needs to be something like maybe 3 for galfit
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']
        gain[i] = header[i]['CCDGAIN']
        RN[i] = header[i]['READNSEA']

#GALFIT RESIDUAL J0905 coarse 475 
        galfile = glob.glob('/Volumes/physics/linux-lab/data/galfit/'+galaxies[2]+'*/J0905_F475W_coarse.fits')
        multi = fits.open(galfile[0])
        multi_res, multi_header = multi[3].data, multi[3].header
        multi_data = multi[1].data
        #define positions for photometry
        positions = [(xcen[w], ycen[w])]
        # do pixel analysis
        for j in range(0,width):
            for k in range(0,width):
                fluxpix[w,i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))
                SNR[w,i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                pixNoise[w,i,width-j-1,k] =  math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i])
m = np.ma.masked_where(SNR<10,SNR)

#YOOOOOO lipscomb radius yo 
#snrdisoverlap=np.zeros([len(galaxies),len(wavelengths), width, width,2])
#radiigroups=np.empty([12,3,2])
#for w in range(0,len(galaxies)):
#    for z in range(0,len(wavelengths)):
#        snrdisoverlap[w][z] = np.dstack((disfrombest[w],SNR[w][z]))
#for w in range(0,len(galaxies)):
#    for z in range(0,len(wavelengths)):
#        for q in range(0,width):
#            for p in range(0,width):
#                if snrdisoverlap[w][z][q][p][0]<9:
#                    np.append(radiigroups[w][z][0], snrdisoverlap[w][z][q][p])


# Define residual J0905 coarse noise and plot it against radius
snrJ0905 = SNR[2]
nsrJ0905 = 1/snrJ0905
fluxpixJ0905 = fluxpix[2]
noiseJ0905 = nsrJ0905*fluxpixJ0905
residual_noise_J0905 = noiseJ0905*(2**0.5)

#make plot in pdf of J0905 residual noise vs radius from centroid
with PdfPages('cl_residual_noise_J0905.pdf') as pdf:   
    plt.figure()

    plt.scatter(disfrombest[2],residual_noise_J0905[0])
    plt.scatter(disfrombest[2],residual_noise_J0905[1])
    plt.scatter(disfrombest[2],residual_noise_J0905[2])
    plt.ylim(0,200)
    plt.xlim(0,40)
    plt.xlabel('Distance from Centroid (pix)')
    plt.ylabel('Residual Image Noise (electrons)')
    plt.title('Noise by Distance for J0905 Residuals')
    a=plt.scatter(disfrombest[2],residual_noise_J0905[0])
    b=plt.scatter(disfrombest[2],residual_noise_J0905[1])
    c=plt.scatter(disfrombest[2],residual_noise_J0905[2])
    plt.legend((a,b,c,),
           ('475 coarse', '814 coarse', '1600 coarse'),
           scatterpoints=1,
           loc='upper right',
           ncol=3,
           fontsize=8)
    pdf.savefig()
    plt.close()

#cut out 81 by 81 slice of 401 by 401 J0905 475 coarse residual flux and noise arrays
J0905ressig = np.zeros([width,width])
J0905resnoise = np.zeros([width,width])
for i in range(width):
    for q in range(width):
        J0905ressig[i][q] = multi_res[161+i][161+q]
        J0905resnoise = residual_noise_J0905[0]
        
#define snr for J0905 475 coarse residual
SNRresJ0905 = J0905ressig/J0905resnoise

#make plot of SNRresJ0905 coarse 475 vs radius
with PdfPages('cl_SNRresJ0905 vs radius.pdf') as pdf:   
    plt.figure()
    for k in range(width):
        for j in range(width):
            plt.scatter(disfrombest[2][j][k],SNRresJ0905[j][k])
    plt.ylim(0,100)
    plt.xlim(0,60)
    plt.xlabel('Distance from Centroid (pix)')
    plt.ylabel('SNR Residual Images')
    plt.title('SNR for the Residual Images from J0905 F475W Coarse vs Radius')
    pdf.savefig()
    plt.close()
    
#make plot of signalresJ0905 coarse 475 vs radius
with PdfPages('cl_signalresJ0905 vs radius.pdf') as pdf:   
    plt.figure()
    for k in range(width):
        for j in range(width):
            plt.scatter(disfrombest[2][j][k],J0905ressig[j][k])
    #plt.ylim(0,1000)
    plt.xlim(0,60)
    plt.xlabel('Distance from Centroid (pix)')
    plt.ylabel('Flux in Residual Images (electrons)')
    plt.yscale('log')
    plt.title('Flux of the Residual Images for J0905 F475W Coarse vs Radius')
    pdf.savefig()
    plt.close()
