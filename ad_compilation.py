# Aleks Diamond-Stanic
# 20161122
#
#the goal of this script is to take our flux values, whether it be in annulus form or total flux form, and calculate a few things:
# 1. luminosity distance using Ned Wright's calculator (use calc. online, then write code for it later)
# 2. measure total flux in galaxy
# 3. compute luminosity using Hogg equation 24
# 4. compute mass-to-light ratio using Bell & de Jong Table 1
# 5. compute mass of galaxy
# 6. for one galaxy, measure flux, color, luminosity, mass-to-light ratio and mass for each annulus
# 7. make plots of flux vs. radius, color vs. radius, and mass vs. radius for a galaxy

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import math

# define the directory that contains the images
dir = os.environ['HSTDIR']

#setting up arrays with three elements, all zeros - placeholders
wavelengths = [4,8,1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3388.
ycen = 3504.
dx = 100
dy = 100

# define the radii to be used for aperture photometry
radii = np.arange(40)+1

#make an array for the calculation of the area of each bagel (annulus)
area = [0 for x in range(len(radii))]

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# create a PDF file for the plots    
with PdfPages('ad_compilation.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir+'final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']

        #define positions for photometry
        positions = [(xcen, ycen)]

        #do photometry on images
        #convert to proper units
        for j in range(0,len(radii)):
            aperture = CircularAperture(positions, radii[j])
            phot_table = aperture_photometry(data[i], aperture)
            flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])
            if j == 0:
                subflux[i,j] = flux[i,j]
            else:
                subflux[i,j] = flux[i,j]-flux[i,j-1]
        
        #set up the plot
        #ax = fig.add_subplot(1,1,1)
        #cax = ax.scatter(-2.5*np.log10(subflux[1] / subflux[2]), -2.5*np.log10(subflux[0] / subflux[1]), c=radii, vmin=radii[0], vmax=radii[-1], cmap=cm.coolwarm, s=25, lw=0.2,)

        
    #calculating galaxy-wide
    
    #finding total flux in galaxy in erg/s
    tflux475 = flux[0,len(radii)-1]*(10**-23)
    tflux814 = flux[1,len(radii)-1]*(10**-23)
    tflux160 = flux[2,len(radii)-1]*(10**-23)

    #luminosity distance (in cm) for J0905, z = 0.712
    LdJ0905 = 4348.9*3.08568*10**24

    #finding magnitudes and color for M/L ratio
    mag475 = -2.5*np.log10(tflux475 / 3631)
    mag814 = -2.5*np.log10(tflux814 / 3631)
    mag160 = -2.5*np.log10(tflux160 / 3631)
    colorUV = mag475-mag814
    colorVJ = mag814-mag160

    #determining M/L ratio using Table 1 of Bell & de Jong
    MLR_BV_B = 10**(-.994+(1.804*colorUV))
    MLR_BV_V = 10**(-.734+(1.404*colorUV))
    MLR_BV_J = 10**(-.621+(.794*colorUV))
    MLR_VJ_B = 10**(-1.903+(1.138*colorVJ))
    MLR_VJ_V = 10**(-1.477+(.905*colorVJ))
    MLR_VJ_J = 10**(-1.029+(.505*colorVJ))

    #calculating nu_e * L_nu_e luminosity in erg/s units from Hogg eq (24)
    c = 299792458
    LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0905**2)
    LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0905**2)
    LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0905**2)

    #convert luminosity to solar units
    Lsol475 = LnuNu475 / (3.846*10**33)
    Lsol814 = LnuNu814 / (3.846*10**33)
    Lsol160 = LnuNu160 / (3.846*10**33)

    #calculate mass of galaxy in solar units
    M_BV_B = Lsol475*MLR_BV_B
    M_BV_V = Lsol814*MLR_BV_V
    M_BV_J = Lsol160*MLR_BV_J
    M_VJ_B = Lsol475*MLR_VJ_B
    M_VJ_V = Lsol814*MLR_VJ_V
    M_VJ_J = Lsol160*MLR_VJ_J
    print('M_BV_B', M_BV_B/1e11, 'M_BV_V', M_BV_V/1e11, 'M_BV_J', M_BV_J/1e11, 'M_VJ_B', M_VJ_B/1e11, 'M_VJ_V', M_VJ_V/1e11, 'M_VJ_J', M_VJ_J/1e11)
    

    #calculation annulus-based

    #calculation of flux for each annulus, in erg/s units
    aflux475 = subflux[0]*10**-23
    aflux814 = subflux[1]*10**-23
    aflux160 = subflux[2]*10**-23

    #calculation of magnitudes and color for each annulus
    amag475 = -2.5*np.log10(aflux475 / 3631)
    amag814 = -2.5*np.log10(aflux814 / 3631)
    amag160 = -2.5*np.log10(aflux160 / 3631)
    acolorUV = amag475-amag814
    acolorVJ = amag814-amag160

    #determining M/L ratio using Table 1 of Bell & de Jong
    aML_BV_B = 10**(-.994+(1.804*acolorUV))
    aML_BV_V = 10**(-.734+(1.404*acolorUV))
    aML_BV_J = 10**(-.621+(.794*acolorUV))

    #calculating nu_e * L_nu_e luminosity in erg/s units for each annulus from Hogg eq (24)
    c = 299792458
    aLnuNu475 = (c/(475*10**-9))*aflux475*(4*math.pi*LdJ0905**2)
    aLnuNu814 = (c/(814*10**-9))*aflux814*(4*math.pi*LdJ0905**2)
    aLnuNu160 = (c/(1600*10**-9))*aflux160*(4*math.pi*LdJ0905**2)

    #convert luminosity for each annulus to solar units
    aLsol475 = aLnuNu475 / (3.846*10**33)
    aLsol814 = aLnuNu814 / (3.846*10**33)
    aLsol160 = aLnuNu160 / (3.846*10**33)

    #calculate mass associated with each annulus in solar units
    aM_BV_B = aLsol475*aML_BV_B
    aM_BV_V = aLsol814*aML_BV_V
    aM_BV_J = aLsol160*aML_BV_J
            

    # set plot parameters
    #cbar = fig.colorbar(cax)
    #cbar.set_label('radius [pixels]', fontsize=18)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #plt.ylim([10**(-1),1e1])
    #plt.ylabel('U-V', fontsize=18)
    #plt.xlim([10**(-1),1e1])
    #plt.xlabel('V-J', fontsize=18)

    pdf.savefig()
    plt.close()

    
    #os.system('open %s &' % 'ad_compilation.pdf')
