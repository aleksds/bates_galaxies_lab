# Sophia C W Gottlieb I
# 20161202
#
# This code shall take jr_compilation and jr_overlay and make them into one.
# Another thing we intend this to do is to run all galaxies in a loop.
# Is this possible? We shall see....

# Almost all of this code was originally written by J. H. Rines in 201611. 
#
# The goals of J H R's code were as follows:
# take our flux values, whether it be in annulus form or total flux form, and calculate a few things:
# 1. luminosity distance using Ned Wright's calculator (use calc. online, then write code for it later)
# 2. measure total flux in galaxy
# 3. compute luminosity using Hogg equation 24
# 4. compute mass-to-light ratio using Bell & de Jong Table 1
# 5. compute mass of galaxy
# 6. for galaxy J0905, measure flux, color, luminosity, mass-to-light ratio and mass for each annulus
# 7. make plots of flux vs. radius, color vs. radius, and mass vs. radius for a galaxy
#
# In addition to those goals previously mentioned, this script shall do this for many galaxies.
#
# Things that my brain is thinking to do.....
# 1. find and replace J0905 with galaxy[0].
# 2. find the discrepancies between comp and overlay to create a successful merge.
# 3. take out all hard coded numbers, replace with descriptive variable names
# 4. run this mother function in a loop.
# 5. Maybe cut this up into functions?
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
from matplotlib.backends.backend_pdf import PdfPages
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from scipy import integrate

# define the directory that contains the images
dir = os.environ['HSTDIR']

#setting up arrays with three elements, all zeros - placeholders
wavelengths = [4,8,1]
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]

galaxies = ['J0905']

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3386.5
ycen = 3503.2
dx = 100
dy = 100

# define the radii to be used for aperture photometry
radii = np.arange(40)+1

# make an array for the calculation of the area of each bagel (annulus)
area = [0 for x in range(len(radii))]

# calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# create a PDF file for the plots
with PdfPages('sg_comp_overlay_' + galaxies[0] +'.pdf') as pdf:

    fig = plt.figure()
 
    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir + galaxies[0] + '_final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']

        # define positions for photometry
        positions = [(xcen, ycen)]

        # do photometry on images
        # convert to proper units
        for j in range(0,len(radii)):
            aperture = CircularAperture(positions, radii[j])
            phot_table = aperture_photometry(data[i], aperture)
            flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])
            if j == 0:
                subflux[i,j] = flux[i,j]
            else:
                subflux[i,j] = flux[i,j]-flux[i,j-1]

    # calculating galaxy-wide
# THIS IS WHERE GOTTLIEB STOPS UNDERSTANDING
# Look up and define new variables
    # set up two-dimensional arrays for the a and b coefficients based on the luminosity and color
    # this will be in the same format ish as the table in Josh's blue notebook
    # Bell 2001 Table 3
    Ba = [-1.019,-1.113,-1.026,-.990,-1.110,-.994,-.888]
    Bb = [1.937,2.065,1.954,1.883,2.018,1.804,1.758]
    B_coeff = [Ba,Bb]
    Va = [-.759,-.853,-.766,-.730,-.850,-.734,-.628]
    Vb = [1.537,1.665,1.554,1.483,1.618,1.404,1.358]
    V_coeff = [Va,Vb]
    Ja = [-.540,-.658,-.527,-.514,-.659,-.621,-.550]
    Jb = [.757,.907,.741,.704,.878,.794,.801]
    J_coeff = [Ja,Jb]
    
    # finding total flux in galaxy in erg/s
    tflux475 = flux[0,len(radii)-1]*(10**-23)
    tflux814 = flux[1,len(radii)-1]*(10**-23)
    tflux160 = flux[2,len(radii)-1]*(10**-23)

    # luminosity distance (in cm) for galaxy, z = 0.712
# This will have to be specified for each galaxy
# I don't understand where all of these numbers came from
# From Ned Wright online Calculator
    LdJ0905 = 4348.9*3.08568*10**24

    # finding magnitudes and color for M/L ratio
# What is this -2.5 and this 3631?
# -2.5 
    mag475 = -2.5*np.log10(tflux475 / 3631)
    mag814 = -2.5*np.log10(tflux814 / 3631)
    mag160 = -2.5*np.log10(tflux160 / 3631)
    colorUV = mag475-mag814
    colorVJ = mag814-mag160

    # determining M/L ratio using Table 3 of Bell & de Jong, seven coefficients for each luminosity for BV color
# Mass to Light ratio, but then I'm not sure what's happening
    MLR_BV_Bk = np.zeros([len(Ba)])
    MLR_BV_Vk = np.zeros([len(Va)])
    MLR_BV_Jk = np.zeros([len(Ja)])
    for k in range(0,len(Ba)):
        MLR_BV_Bk[k] = 10**(Ba[k]+(Bb[k]*colorUV))
        MLR_BV_Vk[k] = 10**(Va[k]+(Vb[k]*colorUV))
        MLR_BV_Jk[k] = 10**(Ja[k]+(Jb[k]*colorUV))

    # calculating nu_e * L_nu_e luminosity in erg/s units from Hogg eq (24), only three values depending on filter
# Calculates *L nu? Lnu is Lumonosity? c is sthe speed of light... speed / frequency *total fluc * surface area?
    c = 299792458
    LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0905**2)
    LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0905**2)
    LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0905**2)

    #convert luminosity to solar units
# 3.846*10**33 is the mass of the sun
    Lsol475 = LnuNu475 / (3.846*10**33)
    Lsol814 = LnuNu814 / (3.846*10**33)
    Lsol160 = LnuNu160 / (3.846*10**33)

    #calculate mass of galaxy in solar units, 3 arrays of masses for each coefficient for 3 different filters will yield a total of 21 galaxy masses, 7 values for each luminosity in the BV color
    M_BV_Bk = np.zeros([len(Ba)])
    M_BV_Vk = np.zeros([len(Va)])
    M_BV_Jk = np.zeros([len(Ja)])
    
    M_BV_Bk = Lsol475*MLR_BV_Bk
    M_BV_Vk = Lsol814*MLR_BV_Vk
    M_BV_Jk = Lsol160*MLR_BV_Jk

    mass = (M_BV_Bk, M_BV_Vk, M_BV_Jk)

    #calculation annulus-based

    #calculation of flux for each annulus, given in an array, for each filter, in erg/s units
    aflux475 = subflux[0]*10**-23
    aflux814 = subflux[1]*10**-23
    aflux160 = subflux[2]*10**-23

    #calculation of magnitudes and color for each annulus
    #I do not think I need to calculate 'acolorUV', but rather simply use the 'colorUV' calculation to determine MLR, however I have left this in the code in case we need to refer to it later. I did not comment it out since it isn't used later in the code.
    amag475 = -2.5*np.log10(aflux475 / 3631)
    amag814 = -2.5*np.log10(aflux814 / 3631)
    amag160 = -2.5*np.log10(aflux160 / 3631)
    acolorUV = amag475-amag814
    acolorVJ = amag814-amag160

    #I THINK I FIGURED IT OUT: I only need the ML ratio for the entire galaxy, NOT for each specific annulus...that somehow seems to mess up the total mass.  THUS: I am essentially using the same MLR as MLR_BV_X , so the 'acolor' things were not necessary I dont think...but if need be, I will put the 'acolorUV' code back in place of the 'colorUV' code below.
    
    #this is determining MLR using Table 1 coeffieicnts.  I could extend this to use all the coefficients, but for now, my plan is to calculate statistical uncertainty using the 21 values I have from finding the total mass of the galaxy and extend it to these measurements.  I could use the 3 values I have calculated below for the annuli.

    #after discussing with Aleks, we determined that we do want an individualized MLR for each annulus, so I have done that by replacing the 'acolorUV' below. This gives an array with size 40 with which we can calculate the aLunuNu, and later on the mass in each ring.
    aMLR_BV_B = 10**(-.994+(1.804*acolorUV))
    aMLR_BV_V = 10**(-.734+(1.404*acolorUV))
    aMLR_BV_J = 10**(-.621+(.794*acolorUV))

    #calculation of annular MLR based on one MLR for the entire galaxy, since we plan to overlay polots of both. this will be denoted 'bMLR...'
    bMLR_BV_B = 10**(-.994+(1.804*colorUV))
    bMLR_BV_V = 10**(-.734+(1.404*colorUV))
    bMLR_BV_J = 10**(-.621+(.794*colorUV))

    #calculating nu_e * L_nu_e luminosity in erg/s units for each annulus from Hogg eq (24)
    c = 299792458
    aLnuNu475 = (c/(475*10**-9))*aflux475*(4*math.pi*LdJ0905**2)
    aLnuNu814 = (c/(814*10**-9))*aflux814*(4*math.pi*LdJ0905**2)
    aLnuNu160 = (c/(1600*10**-9))*aflux160*(4*math.pi*LdJ0905**2)

    #NOTE: we can use the same luminosity calculations, since the luminosity doesnt depend on M/L ratio, which is why we don't make a 'bLunuNu...' calculation
    
    #convert luminosity for each annulus to solar units
    aLsol475 = aLnuNu475 / (3.846*10**33)
    aLsol814 = aLnuNu814 / (3.846*10**33)
    aLsol160 = aLnuNu160 / (3.846*10**33)

    #calculate mass associated with each annulus, based on individualized MLRs for each annulus, in solar units, based on individualized MLRs
    aM_BV_B = aLsol475*aMLR_BV_B
    aM_BV_V = aLsol814*aMLR_BV_V
    aM_BV_J = aLsol160*aMLR_BV_J
    amass = (aM_BV_B,aM_BV_V,aM_BV_J)

    #calculate mass associated with each annulus in solar units, based on one MLR estimate for the entire galaxy
    bM_BV_B = aLsol475*bMLR_BV_B
    bM_BV_V = aLsol814*bMLR_BV_V
    bM_BV_J = aLsol160*bMLR_BV_J
    bmass = (bM_BV_B,bM_BV_V,bM_BV_J) 

    #making various calculations, including std, best value, etc.
    ave_mass475 = np.mean(mass[0])
    ave_mass814 = np.mean(mass[1])
    ave_mass160 = np.mean(mass[2])
    ave_mass = (ave_mass475,ave_mass814,ave_mass160)
    
    std_mass475 = np.std(mass[0])
    std_mass814 = np.std(mass[1])
    std_mass160 = np.std(mass[2])
    std_mass = (std_mass475,std_mass814,std_mass160)

    step1 = ave_mass[0]/np.square(std_mass[0]) + ave_mass[1]/np.square(std_mass[1]) + ave_mass[2]/np.square(std_mass[2])
    step2 = 1/np.square(std_mass[0]) + 1/np.square(std_mass[1]) + 1/np.square(std_mass[2])

    best_value = step1/step2
    print('J0905 total mass best value:', best_value/1e11)

    uncert = (1/np.square(std_mass[0]) + 1/np.square(std_mass[1]) + 1/np.square(std_mass[2]))**-.5
    print('J0905 total mass uncertainty:', uncert/1e11)

    #plotting mass vs radius
    acolors = ['b--','g--','r--']
    bcolors = ['b', 'g', 'r']
    adot = ['b','g','r']
    bdot = ['bo','go','ro']
    alabeling = ['annular MLR F475W','annular MLR F814W','annular MLR F160W']
    blabeling = ['single MLR F475W','single MLR F814W','single MLR F160W']


    # THIS IS WHERE OVERLAY AND COMP DIVERGE. Line 269 in jr_comp
    
    #putting radius into kpc
    kpc_radius = radii*(0.05)*(7.194)
# COMP has something similar, but slightly different 
    #plotting the specific annular (specific aMLR) mass
    ax = fig.add_subplot(2,1,1)
    ax.plot(kpc_radius, amass[1], acolors[1], marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, amass[1], acolors[1])
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.title('J0905 Mass vs. Radius, annular M/L ratios',fontsize=16)
    plt.tight_layout()
    legend = ax.legend(loc='upper right')
        
    #plotting the broad (single bMLR) annular mass, the 'light mass' as I call it
    
    bx = fig.add_subplot(2,1,1)
    bx.plot(kpc_radius, bmass[1], 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bmass[1], 'yellowgreen')

    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.tight_layout()
    plt.title(galaxies[0] + ' Mass vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')

    
    #plotting the mass/area vs radius, converting radius units from pixels to kpc to get mass surface density units
    
    #first creating an array with areas of shells in proper units of kpc
    kpc_area = np.zeros(len(area))
    amass1_ovr_area = np.zeros(len(area))
    bmass1_ovr_area = np.zeros(len(area))
    for j in range(0,len(area)):
       kpc_area[j] = area[j]*(0.05**2)*(7.194**2)
       
    #now calculating mass/area (mass surface density) in units of solar masses/kpc for amass and bmass and also radius in kpc units
    amass1_ovr_area = amass[1]/kpc_area
    bmass1_ovr_area = bmass[1]/kpc_area
    kpc_radius = radii*(0.05)*(7.194)
    
    #now plotting amass1_ovr_area and bmass1_ovr_area vs radius in kpc
    ax = fig.add_subplot(2,1,2)
    ax.plot(kpc_radius, amass1_ovr_area, 'g--', marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, amass1_ovr_area, 'g--')
   
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title('J0905 Mass Density vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')

    #plotting bmass1_ovr_area vs radius in kpc
    bx = fig.add_subplot(2,1,2)
    bx.plot(kpc_radius, bmass1_ovr_area, 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bmass1_ovr_area, 'yellowgreen')
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title(galaxies[0] + ' Mass Density vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')
    
    #calculating total mass (amass) for annular MLR (814 filter only)
    total_annular_amass_F814W = np.sum(amass[1])
    print('total amass', total_annular_amass_F814W)

    #calculating total mass (bmass) for single MLR (814 filter only)
    total_singular_bmass_F814W = np.sum(bmass[1])
    print('total bmass', total_singular_bmass_F814W)


    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_overlays_' + galaxies[0] + '.pdf')
