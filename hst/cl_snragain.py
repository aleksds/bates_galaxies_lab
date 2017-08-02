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

# define the directory that contains the images
dir = os.environ['HSTDIR']
#dir = 'aaw/data/'

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


# Now, we loop through all galaxies
with PdfPages('cl_snrcoarsecolor.pdf') as pdf:
    for w in range (0, len(galaxies)):
                
        for i in range (0, 3):
            
            # read in the images
            file = glob.glob(dir+galaxies[w]+'*/coarse/'+wavelengths[i]+'/final*sci.fits')
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
                    fluxpix[w,i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))
                    SNR[w,i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                    pixNoise[w,i,width-j-1,k] =  math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i])
   
        acolors = ['b--','g--','r--','c--']
        bcolors = ['b', 'g', 'r']
        dot = ['bo','go','ro']
        labeling = ['475 nm','814 nm','1600 nm']
        siglabel = ['read','sky','dark','source']

        m = np.ma.masked_where(SNR<10,SNR)
        for i in range(0,3):
            fig = plt.figure()
            bx = fig.add_subplot(1,1,1)
            plt.imshow(m[w][i])
            plt.colorbar()
            plt.title(galaxies[w]+ ' SNR at ' +str(wavelengths[i]))
            plt.xlabel('Pixels')
            plt.ylabel('Pixels')
            pdf.savefig()
            plt.close()

#from cl_centroidmerged.py

bestxswhole = np.array([ 3628.8630965 ,  3933.65237012,  3386.54069503,  3477.20367975, 3572.97219824,  3803.17972872,  3883.20568985,  4147.00201888, 3786.77909178,  4175.294108  ,  3565.96733481,  4065.88540274])
bestyswhole = np.array([ 4153.8967805 ,  4136.5722187 ,  3503.1663102 ,  3405.39675221, 3339.1874928 ,  4170.32537204,  4164.77089952,  3922.116162  , 4186.05759774,  3826.62884958,  3434.92868183,  4053.58473705])

xcendiff = bestxswhole-xcen
ycendiff = bestyswhole-ycen

disfrombest=np.zeros([len(galaxies),(width),(width)])
punc = (pixNoise/fluxpix)*100
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
                disfrombest[w][p-1][q-1] = ((41-np.abs(xcendiff[w]))**2 + (41-q+np.abs(ycendiff[w]))**2)**0.5
            if p > 41 and q == 41:
                disfrombest[w][p-1][q-1] = ((41-p-xcendiff[w])**2 + (41-q+np.abs(ycendiff[w]))**2)**0.5
            if p < 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p+xcendiff[w])**2 + (41-q+ycendiff[w])**2)**0.5
            if p == 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p+np.abs(xcendiff[w]))**2 + (41-q+ycendiff[w])**2)**0.5
            if p > 41 and q > 41:
                disfrombest[w][p-1][q-1] = ((41-p-xcendiff[w])**2 + (41-q+ycendiff[w])**2)**0.5
    #disfrombest = disfrombest*pixtokpc[w]
with PdfPages('cl_snrcourserad.pdf') as pdf:
    for w in range(0,len(galaxies)):
        fig = plt.figure()
        plt.subplots_adjust(wspace=.6,hspace=.6)
        ax = fig.add_subplot(2,2,1)
        bx = fig.add_subplot(2,2,2)
        cx = fig.add_subplot(2,2,3)
        ax.scatter(disfrombest[w],SNR[w][0])
        bx.scatter(disfrombest[w],SNR[w][1])
        cx.scatter(disfrombest[w],SNR[w][2])

        ax.set_title(galaxies[w]+' F475W',fontsize=10)
        ax.set_ylabel('SNR')
        ax.set_xlabel('radius (pix)')
        ax.set_ylim(0,20)
        bx.set_title(galaxies[w]+' F814W',fontsize=10)
        bx.set_ylabel('SNR')
        bx.set_xlabel('radius (pix)')
        bx.set_ylim(0,20)
        cx.set_title(galaxies[w]+' F160W',fontsize=10)
        cx.set_ylabel('SNR')
        cx.set_xlabel('radius (pix)')
        cx.set_ylim(0,20)
        pdf.savefig()
        plt.subplots_adjust(wspace=10,hspace=10)
        plt.close()
        
#FINEDATA
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J2116', 'J2140']
xcen = [7260,7870,6775,6956,7148,7609,7769,8296,7576,7134,8134]
ycen = [8310,8276,7009,6813,6680,8343,8332,7846,8374,6872,8110]



finewidth = 161
pixr = 80
fluxpix = np.zeros([len(galaxies),len(wavelengths), finewidth, finewidth])
pixNoise = np.zeros([len(galaxies),len(wavelengths), finewidth, finewidth])
SNR = np.zeros([len(galaxies),len(wavelengths), finewidth, finewidth])

with PdfPages('cl_snrfinecolor.pdf') as pdf:
    for w in range (0, len(galaxies)):
                
        for i in range (0, 2):
            
            # read in the images
            file = glob.glob(dir+galaxies[w]+'*/fine/'+wavelengths[i]+'/final*sci.fits')
            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
            RN[i] = header[i]['READNSEA']

            #define positions for photometry
            positions = [(xcen[w], ycen[w])]
            # do pixel analysis
            for j in range(0,finewidth):
                for k in range(0,finewidth):
                    fluxpix[w,i,finewidth-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))
                    SNR[w,i,finewidth-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                    pixNoise[w,i,finewidth-j-1,k] =  math.sqrt((rSky[i][w]*1)**2+fluxpix[w][i][finewidth-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i])
   
        acolors = ['b--','g--','r--','c--']
        bcolors = ['b', 'g', 'r']
        dot = ['bo','go','ro']
        labeling = ['475 nm','814 nm','1600 nm']
        siglabel = ['read','sky','dark','source']

        m = np.ma.masked_where(SNR<10,SNR)
        for i in range(0,2):
            fig = plt.figure()
            bx = fig.add_subplot(1,1,1)
            plt.imshow(m[w][i])
            plt.colorbar()
            plt.title(galaxies[w]+ ' SNR at ' +str(wavelengths[i]))
            plt.xlabel('Pixels')
            plt.ylabel('Pixels')
            pdf.savefig()
            plt.close()

#from cl_centroidforfine.py
bestxswhole = np.array([ 7258.63198319,  7868.47231462,  6774.20583264,  6955.46121619, 7146.80768432,  7607.48957592,  7767.61187256,  8294.97520887, 7574.55587795,  7132.91470471,  8132.85741776])
bestyswhole = np.array([ 8308.96520047,  8274.48208819,  7007.36845459,  6811.81892727, 6679.3589637 ,  8341.98082999,  8330.66769695,  7845.37390573, 8373.24262764,  6870.9730902 ,  8108.24043305])

xcendiff = bestxswhole-xcen
ycendiff = bestyswhole-ycen

disfrombest=np.zeros([len(galaxies),(finewidth),(finewidth)])

#DIFFERENT FOR FINE
for w in range(0,11):
    for q in range(1,finewidth+1):
        for p in range(1,finewidth+1):
            if p < 81 and q < 81:
                disfrombest[w][p-1][q-1] = ((81-p+xcendiff[w])**2 + (81-q-ycendiff[w])**2)**0.5
            if p == 81 and q < 81:
                disfrombest[w][p-1][q-1] = ((81-p+np.abs(xcendiff[w]))**2 + (81-q-ycendiff[w])**2)**0.5
            if p > 81 and q < 81:
                disfrombest[w][p-1][q-1] = ((81-p-xcendiff[w])**2 + (81-q-ycendiff[w])**2)**0.5
            if p < 81 and q == 81:
                disfrombest[w][p-1][q-1] = ((81-p+xcendiff[w])**2 + (81-q+np.abs(ycendiff[w]))**2)**0.5
            if p == 81 and q == 81:
                disfrombest[w][p-1][q-1] = ((81-np.abs(xcendiff[w]))**2 + (81-q+np.abs(ycendiff[w]))**2)**0.5
            if p > 81 and q == 81:
                disfrombest[w][p-1][q-1] = ((81-p-xcendiff[w])**2 + (81-q+np.abs(ycendiff[w]))**2)**0.5
            if p < 81 and q > 81:
                disfrombest[w][p-1][q-1] = ((81-p+xcendiff[w])**2 + (81-q+ycendiff[w])**2)**0.5
            if p == 81 and q > 81:
                disfrombest[w][p-1][q-1] = ((81-p+np.abs(xcendiff[w]))**2 + (81-q+ycendiff[w])**2)**0.5
            if p > 81 and q > 81:
                disfrombest[w][p-1][q-1] = ((81-p-xcendiff[w])**2 + (81-q+ycendiff[w])**2)**0.5

with PdfPages('cl_snrfinerad.pdf') as pdf:
    for w in range(0,len(galaxies)):
        fig = plt.figure()
        plt.subplots_adjust(wspace=.6,hspace=.6)
        ax = fig.add_subplot(2,2,1)
        bx = fig.add_subplot(2,2,2)

        ax.scatter(disfrombest[w],SNR[w][0])
        bx.scatter(disfrombest[w],SNR[w][1])
        
        ax.set_title(galaxies[w]+' F475W',fontsize=10)
        ax.set_ylabel('SNR')
        ax.set_xlabel('radius (pix)')
        ax.set_ylim(0,20)
        ax.set_xlim(0,80)
        bx.set_title(galaxies[w]+' F814W',fontsize=10)
        bx.set_ylabel('SNR')
        bx.set_xlabel('radius (pix)')
        bx.set_ylim(0,20)
        bx.set_xlim(0,80)
        
        pdf.savefig()
        plt.subplots_adjust(wspace=10,hspace=10)
        plt.close()
