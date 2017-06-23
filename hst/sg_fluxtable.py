# Sophia C W Gottlieb I
# 20170612 produce table of flux values in a txt file
# This code was originally the sg_compoverlay_loop.py, but I am going to take all of that code out right about now.
#
# import relevant Python modules
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

# define the directory that contains the images
dir = os.environ['HSTDIR']

#setting up arrays with three elements, all zeros - placeholders
wavelengths = [4,8,1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]
gain = [0 for x in range(len(wavelengths))]

#width MUST be odd.
width = 15
pixr = int((width-1)/2)

# specify the position of the science target and the size of the region around the science target to consider
filters = np.array([475, 814, 1600]) #*u.nm
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
xcen = [3628, 3933, 3386.5, 3477.5, 3573, 3802, 3886, 4149, 3787, 4174, 3565, 4067]
ycen = [4153, 4136, 3503.2, 3404.3, 3339, 4169, 4164, 3921, 4187, 3826, 3434, 4054]

fluxvalues = [[0 for x in range(len(wavelengths))] for y in range(len(galaxies))]

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
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# Now, we loop through all galaxies
for w in range (0, len(galaxies)):
#for w in range (2,3):
    # Mostly just for our own benefit, we print the galaxy's name so we know what 
    # works and what doesn't when the code breaks.
    print(galaxies[w])
    with PdfPages('sg_flux_'+galaxies[w]+'.pdf') as pdf: 
        fig = plt.figure()
        collection = ['F475W','F814W','F160W']
    
        flux = np.zeros([len(collection),len(radii)]) #*u.Jy
        subflux = np.zeros([len(collection),len(radii)])
        fluxpix = np.zeros([len(collection), width, width])
        for i in range (0, len(collection)):
            
            # read in the images
            file = glob.glob(dir+galaxies[w]+'_final_'+collection[i]+'*sci.fits')
            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
    
            #define positions for photometry
            positions = [(xcen[w], ycen[w])]
    
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    #print(i,j,k)
                    #fluxpix[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*gain[i]*exp[i])
                    fluxpix[i,width-j-1,k] = math.log10((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*gain[i]*exp[i])
                    totalphotons = totalphotons +  (data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*gain[i]*exp[i]
                    #converts units to Jy and then to nanomaggys: Jy is data * fnu / exp and 1 nm = 3.631e-6 Jy
                    #fluxnmaggys[i,width-j-1,k] = data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*fnu[i]/exp[i]/(3.631*10**-6)

            plt.imshow(fluxpix[i],cmap='gray')
            plt.colorbar()
            plt.title(galaxies[w]+' photon count: '+str(totalphotons) +' ('+ collection[i] + ')')
            plt.xlabel('Pixels')
            plt.ylabel('Pixels')
            pdf.savefig()
            plt.close()
# NANOMAGGYS
            for j in range(0,len(radii)):
                aperture = CircularAperture(positions, radii[j])
                phot_table = aperture_photometry(data[i], aperture)
                # CONVERT UNITS TO JY TO nMy
                flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])/(3.631*10**-6)
                #flux[i,j] = phot_table['aperture_sum'][0]*gain[i]*exp[i]
                if j == 0:
                    subflux[i,j] = flux[i,j]
                else:
                    subflux[i,j] = flux[i,j]-flux[i,j-1]

            fluxvalues[w]=subflux

# Now that fluxvalues is full, we have to start making a new table.....
# the table has 8 columns: ID, f475, iv475, f814, iv814, f1600, iv1600, and z
# ID is GalaxyName_Aperature and all data is in nanomaggys

# we create or open our txt file
f = open("sg_fluxtable_nm.txt","w+")

# we write our column titles - unsure if these need to stay, just thought it would be nice.
f.write('ID\t\tf_475\t\tivar_475\t\tf_814\t\tivar_814\t\tf_1600\t\tivar_1600\t\tz\n')

# we shift through the data by galaxy and then by aperture.
for w in range(0,len(galaxies)):
    for i in range(0,len(radii)):
        # Building the ID name using if/else for those with single digits.
        ID = galaxies[w]+'_'
        if i < 10:
            ID = ID + '0' + str(i)
        else:
            ID = ID + str(i)
        # writing data divided by tabs (\t character)
        f.write(ID+'\t')
        # we loop through the filters because i am lazy and it is technically good form.
        # at the same time, we grab the inverse variance (squared) which is 1/(flux*res)^2.
        for j in range(0,len(collection)):
            f.write(str(fluxvalues[w][j][i])+'\t')
            ivar = (fluxvalues[w][j][i]*percUnc)**(-2)
            f.write(str(ivar)+'\t')
        # to conclude the line, we include the z value for the galaxy before calling for a new line (\n).
        f.write(str(zs[w])+'\n')

f.close()
