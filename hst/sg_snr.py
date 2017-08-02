#Sophia C W Gottlieb I
#June 23 2017
#
# SNR code... The objective is to plot flux v radius with
# error bars using SNR from the PDF in the summer folder
# it also analyzes this on a pixel to pixel basis and applies
# a mask for SNR < 10 to be killed. 
# I would like to apply the mask to the MLR code... more to 
# come on that later.... 

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
xcen = [3628, 3933, int(3386.5), 3477, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycen = [4153, 4136, int(3503.2), 3404, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]
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
fluxpix = np.zeros([len(wavelengths), width, width])

#fluxvalues = [[0 for x in range(len(wavelengths))] for y in range(len(galaxies))]
pixNoise = np.zeros([len(wavelengths), width, width])
SNR = np.zeros([len(wavelengths), width, width])
annNoise = np.zeros([len(wavelengths),len(radii)])

#Ldcm is the luminosity distance in cm, even though astropy thinks it is in Mpc. 
Ldcm = cosmo.luminosity_distance(zs)*u.Mpc.to(u.cm) / u.Mpc
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
radToKpc = conv.arcsec_per_kpc_proper(zs)*0.05/u.arcsec*u.kpc
totalphotons = 0

# define the radii to be used for aperture photometry
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]

# percunc specifies the percent of variation we expect from systematic error... 
# For now, we have chosen 0.05, or 5%.
percUnc = 0.05

#calculate area of each bagel
for i in range(0, len(area)):
    area[i] = math.pi*math.pow(radii[i],2)
    #if i == 0:
    #    area[i] = math.pi*math.pow(radii[0],2)
    #else:
    #    area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# Now, we loop through all galaxies
#for w in range (0, len(galaxies)):
with PdfPages('sg_SNR_err.pdf') as pdf:
    #for w in range(8,9):
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

            #define positions for photometry
            positions = [(xcen[w], ycen[w])]
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    fluxpix[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))#*gain[i]*exp[i])
                    #totalphotons = totalphotons +  (data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])#*gain[i]*exp[i]
                    SNR[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))/(math.sqrt((rSky[i][w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]))
                    pixNoise[i,width-j-1,k] =  math.sqrt((rSky[i][w]*1)**2+fluxpix[i][width-j-1][k]+(RN[i]**2+(gain[i]/2)**2)*1+dark[i]*1*exp[i]) 
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
        
                annNoise[i,j] = math.sqrt((rSky[i][w]*area[j])**2+flux[i][j]+(RN[i]**2)*area[j]+dark[i]*area[j]*exp[i])area = [0 for x in range(len(radii))]
sigskyy = np.zeros([len(wavelengths),len(radii)])
sigread = np.zeros([len(wavelengths),len(radii)])
sigdark = np.zeros([len(wavelengths),len(radii)])
sigstar = np.zeros([len(wavelengths),len(radii)])
                sigskyy[i,j] = (rSky[i][w]*area[j])
                sigread[i,j] = math.sqrt(RN[i]**2+(gain[i]/2)**2*area[j])
                sigdark[i,j] = math.sqrt(dark[i]*area[j]*exp[i])
                sigstar[i,j] = math.sqrt(flux[i][j])
        acolors = ['b--','g--','r--','c--']
        bcolors = ['b', 'g', 'r']
        dot = ['bo','go','ro']
        labeling = ['475 nm','814 nm','1600 nm']
        siglabel = ['read','sky','dark','source']
        #plt.suptitle('Percent Uncertainty for ' +galaxies[w], fontsize=10)
        ax = fig.add_subplot(2,2,1)
        ax.set_title('Total Percent Uncertainty by Filter',fontsize=10)
        
        #bx = fig.add_subplot(2,2,2)
        #cx = fig.add_subplot(2,2,3)
        #dx = fig.add_subplot(2,2,4)
        #ax.plot([radii,radii,radii], annNoise/flux*100,color = 'b', marker='o', label=labeling, linestyle = 'None')
        for k in range(0,len(filters)):
            bx = fig.add_subplot(2,2,k+2)
            bx.set_title('Percent Uncertainty for '+labeling[k]+' by source',fontsize=10)
            ax.plot(radii,(annNoise[k]/flux[k]*100), acolors[k], marker='o', label=str(labeling[k]))
            bx.plot(radii,(sigread[k]/flux[k]*100), acolors[0],marker = 'x',linestyle = 'None',label='read')
            bx.plot(radii,(sigskyy[k]/flux[k]*100), acolors[1],marker = 'o',linestyle = 'None',label='sky')
            bx.plot(radii,(sigdark[k]/flux[k]*100), acolors[2],marker = '^',linestyle = 'None',label='dark')
            bx.plot(radii,(sigstar[k]/flux[k]*100), acolors[3],marker = 's',linestyle = 'None',label='source')
            #ax.plot(radii,np.log10(annNoise[k]/flux[k]*100), acolors[k], marker='o', label=str(labeling[k]))
            #bx.plot(radii,np.log10( sigread[k]/flux[k]*100), acolors[0], marker = 'x',linestyle = 'None',label='read')
            #bx.plot(radii,np.log10( sigskyy[k]/flux[k]*100), acolors[1], marker = 'o',linestyle = 'None',label='sky')
            #bx.plot(radii,np.log10( sigdark[k]/flux[k]*100), acolors[2], marker = '^',linestyle = 'None',label='dark')
            #bx.plot(radii,np.log10( sigstar[k]/flux[k]*100), acolors[3], marker = 's',linestyle = 'None',label='source')
            #plt.ylim(-4,1)
            legend = bx.legend(loc='upper right', fontsize=8)
        
        #bx.plot([radii,radii,radii],sigread/flux*100,color = 'c', linestyle = 'None', marker = 'x', label = 'read')
        #legend = bx.legend(loc = 'upper right')
        plt.xlabel('Radius (pixels)',fontsize=8)
        
        plt.ylabel('Percent Uncertainty',fontsize=8)
        #plt.title(galaxies[w] + ' Percent Uncertainty vs. Radius',fontsize=16)
        plt.tight_layout()
        #legend = plt.legend(loc='upper right')
        legend = ax.legend(loc='upper left')
        # here is some new stuff for you
        pdf.savefig()
        plt.close()   
        #fig = plt.figure()
        m = np.ma.masked_where(SNR<10,SNR)
        for i in range(0,len(filters)):
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
        
# so I think I could make a call like:
# n = np.ma.masked_where(SNR<10, Mass)
# and it would grab all the mass values that work!!!! 
#so let's run this into the other code, and see what falls out!!!!!!!
