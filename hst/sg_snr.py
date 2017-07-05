#Sophia C W Gottlieb I
#June 23 2017
#
# SNR code... The objective is to plot flux v radius with
# error bars using SNR from the PDF in the summer folder

#import relavent packages
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

#width MUST be odd.
width = 81
pixr = int((width-1)/2)
rSky = [ 7.50788545,  7.58130161,  8.27527023,  9.03878993,  8.67763722, 7.62254201,  7.70920672,  6.74143006,  6.87375846,  7.46983987, 7.83102976,  8.10811507]

# specify the position of the science target and the size of the region around the science target to consider
filters = np.array([475, 814, 1600]) #*u.nm
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
xcen = [3628, 3933, int(3386.5), 3477, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycen = [4153, 4136, int(3503.2), 3404, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]

# define the radii to be used for aperture photometry
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]

flux = np.zeros([len(wavelengths),len(radii)]) #*u.Jy
subflux = np.zeros([len(wavelengths),len(radii)])
fluxpix = np.zeros([len(wavelengths), width, width])

#fluxvalues = [[0 for x in range(len(wavelengths))] for y in range(len(galaxies))]
pixNoise = np.zeros([len(wavelengths), width, width])
SNR = np.zeros([len(wavelengths), width, width])
annNoise = np.zeros([len(wavelengths),len(radii)])
zs = [0.603, 0.459, 0.712, 0.514, 0.467, 0.451, 0.658, 0.608, 0.402, 0.449, 0.728, 0.752]

#Ldcm is the luminosity distance in cm, even though astropy thinks it is in Mpc. 
Ldcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc
totalphotons = 0

# define the radii to be used for aperture photometry
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]

# percunc specifies the percent of variation we expect from systematic error... 
# For now, we have chosen 0.05, or 5%.
percUnc = 0.05

#calculate area of each bagel
for i in range(0, len(area)):
    area[i] = math.pi*math.pow(radii[0],2)
    #if i == 0:
    #    area[i] = math.pi*math.pow(radii[0],2)
    #else:
    #    area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# Now, we loop through all galaxies
#for w in range (0, len(galaxies)):
with PdfPages('sg_SNR_err.pdf') as pdf:
    for w in range (0, len(galaxies)):
        print(galaxies[w])
   # with PdfPages('sg_SNR_'+galaxies[w]+'.pdf') as pdf: 
        fig = plt.figure()
                
        for i in range (0, len(wavelengths)):
            
            # read in the images
            file = glob.glob(dir+galaxies[w]+'_final_'+ wavelengths[i]+'*sci.fits')
            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
            RN[i] = header[i]['READNSEA']
            #dark[i] = header[i]['DARKTIME']

            #define positions for photometry
            positions = [(xcen[w], ycen[w])]
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    fluxpix[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))#*gain[i]*exp[i])
                    #totalphotons = totalphotons +  (data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])#*gain[i]*exp[i]
                    SNR[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                    pixNoise[i,width-j-1,k] =  math.sqrt((rSky[w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]) 
                    #converts units to Jy and then to nanomaggys: Jy is data * fnu / exp and 1 nm = 3.631e-6 Jy
                    #fluxnmaggys[i,width-j-1,k] = data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*fnu[i]/exp[i]/(3.631*10**-6)    
            #do photometry on images
            #convert to proper units
            for j in range(0,len(radii)):
                aperture = CircularAperture(positions, radii[j])
                phot_table = aperture_photometry(data[i], aperture)
                #flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])/(3.631*10**-6)
                flux[i,j] = phot_table['aperture_sum'][0]#*gain[i]*exp[i]
                if j == 0:
                    subflux[i,j] = flux[i,j]
                else:
                    subflux[i,j] = flux[i,j]-flux[i,j-1]
                annNoise[i,j] = math.sqrt((rSky[w]*area[j])**2+flux[i][j]+(RN[i]**2+(gain[i]/2)**2)*area[j]+dark[i]*area[j]*exp[i])
        acolors = ['b--','g--','r--']
        bcolors = ['b', 'g', 'r']
        dot = ['bo','go','ro']
        labeling = ['475 nm','814 nm','1600 nm']
        ax = fig.add_subplot(1,1,1)
        for k in range(0,len(acolors)):
            #ax.plot(annNoise[k], flux[k], acolors[k], marker='o', label=str(labeling[k]))
            #ax.errorbar(radii, annNoise[k]/flux[k], yerr = annNoise[k], fmt = 'p')
            ax.plot(radii, annNoise[k]/flux[k]*100, acolors[k], marker='o', label=str(labeling[k]))
            #ax.plot(radii, flux[k], dot[k])
        plt.xlabel('Radius (pixels)',fontsize=14)
        
        plt.ylabel('Percent Uncertainty',fontsize=14)
        plt.title(galaxies[w] + ' Percent Uncertainty vs. Radius',fontsize=16)
        plt.tight_layout()
        #legend = plt.legend(loc='upper right')
        legend = ax.legend(loc='upper right')
        # here is some new stuff for you
        pdf.savefig()
        plt.close()   
        #fig = plt.figure()
        m = np.ma.masked_where(SNR<10,SNR)
        for i in range(0,len(acolors)):
            fig = plt.figure()
            bx = fig.add_subplot(1,1,1)
            plt.imshow(m[i])
            #plt.imshow(SNR[i], cmap='gray')
            #plt.imshow((fluxpix[i]/pixNoise[i]),vmin = 1, vmax = 10,cmap='gray')
            plt.colorbar()
            #plt.imshow(pixNoise[i]/fluxpix[i])
            plt.title(galaxies[w]+ ' SNR at ' +str(wavelengths[i]))
            plt.xlabel('Pixels')
            plt.ylabel('Pixels')
            pdf.savefig()
            plt.close()
        
