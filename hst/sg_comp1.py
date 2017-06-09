# Sophia C W Gottlieb I
# 20170125
#
# This code is the general, optimal version of jr_compilation
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

# set up variables for use
c = 299792458               #c is the speed of light in m/s
solarMass = 3.846*10**33    #solar Mass is the mass of the sun in kg
radToKpc = 0.05*7.194       #converts radius to kpc

# set up experiment specific things
filters = [475, 814, 1600]
galaxies = ['J0905']
xcen = [3386.5]
ycen = [3503.2]
# set up data arrays to be used....
aMass, bMass, aMLR, bMLR = [[0 for x in range(len(filters))] for y in range(len(aB))]
# set up table 1 from Bell & deJong
BellT1 = [-.994, 1.804, -.734, 1.404, -.621, 0.794]
# set up table 3 from Bell & de Jong for Mass to Light Ratio (MLR)
# stellar MLR as a function of color for the scaled salpeter IMF
Ba = [-1.019,-1.113,-1.026,-.990,-1.110,-.994,-.888]
Bb = [1.937,2.065,1.954,1.883,2.018,1.804,1.758]
B_coeff = [Ba,Bb]
Va = [-.759,-.853,-.766,-.730,-.850,-.734,-.628]
Vb = [1.537,1.665,1.554,1.483,1.618,1.404,1.358]
V_coeff = [Va,Vb]
Ja = [-.540,-.658,-.527,-.514,-.659,-.621,-.550]
Jb = [.757,.907,.741,.704,.878,.794,.801]
J_coeff = [Ja,Jb]

# define radii for aperture photometry
radii = np.arange(40)+1

# set up an array for calculation of the area of each bagel (annulus)
# do KPC stuf for later on.
for i in range(0, len(area)):
    kpc_area[i] = area[i]*(radToKpc**2)
    kpc_radius = radii*radToKpc
    
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))


# area calculation for each annulus

# create a PDF file for the plots

with PdfPages('sg_comp1_' + galaxies[0]+ '.pdf') as pdf:
    fig = plt.figure()

    # loop through and collect data from images
    for i in range (0, len(filters)):
        # read in images
        file = glob.glob(dir + galaxies[0] + '_final_'+collection[i]+'*sci.fits')
        hdu = fits.open(file[0])
        data[i], header[i] = hdu[0].data, hdu[0].header
        fnu[i] = header[i]['PHOTFNU']
        exp[i] = header[i]['EXPTIME']
        
        # define photometry positions
        positions = [(xcen[0], ycen[0])]

        # perform photometry on images
        for j in range (0,len(radii)):
            # watch for post-fence error & converts to Jys
            aperture = CircularAperture(positions, radii[j])
            phot_table = aperture_photometry(data[i], aperture)
            flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])
            if j == 0:
                subflux[i,j] = flux[i,j]
            else:
                subflux[i,j] = flux[i,j]-flux[i,j-1]

    # now that the radii and photometry has been completed, we can run some calculations for luminosity
    # specified by wavelength
    for i in range (0, len(filters)):
        # set up each tflux value (total?)
        tflux[i] = flux[i,len(radii)-1]*(10**-23)
        # use tflux for mag
        mag[i] = -2.5*np.log10(tflux[i] / 3631)
        #use tflux for Lsol
        Lsol[i] = (c/(filters[i]*10**-9))*tflux[i]*(4*math.pi*LdJ0905**2) / (solarMass)
        # for every filter, we get an "aflux"? and "amag" where a in ANNULAR.
        aflux[i] = subflux[i]*10**-23
        amag[i] = -2.5*np.log10(aflux[i] / 3631)

    # after that loop, we can now set up color UV by subtracting the magnitude of 475 from 814
    # at the same time, we could set up VJ with 814 from 1600, but this hasn't been used yet
    colorUV = mag[1]-mag[0]
    aColorUV = amag[1]-amag[0]
    # using colorUV, we can determine MLR
    # we loop through all filters, and then the 7 different models used in B&dJ Table 3
    for i in range(0,len(filters)):
        #did josh actually use this? i'm now slightly confused. will print his code and see.
        for j in range (0, len(Ba)):
            # here we calculate MLRs for the 7 models in 3 filters, annularly (aMLR, and broadly bMLR)
            aMLR[i,j] = 10**(aB[j]+(bB[j]*aColorUV))
            # In josh's code, aMLR and bMLR is made specifically with j = 5,
            # which corresponds with Formation epoch: bursts
            # also, josh's code only used i = 1..... meh.
            bMLR[i,j] = 10**(aB[j]+(bB[j]*colorUV))

            # we then get the mass (aMass, bMass) by multiplying the MLR by the L.
            aMass[i,j] = Lsol[i]*aMLR[i,j]
            bMass[i,j] = Lsol[i]*bMLR[i,j]
         # ACTUALLY I DID THIS ELSEWHERE#for every filter, we get an "aflux"? and "amag" where a in ANNULAR.

    # At this point, josh sets up two steps for best value and uncertainty, but i'm not sure he used them
#josh never used them, so i'm not going to set them up right now.
    # we plot mass v radius (CAN I MAKE THIS BETTER?)
    acolors = ['b--','g--','r--']
    bcolors = ['b', 'g', 'r']
    adot = ['b','g','r']
    bdot = ['bo','go','ro']
    alabeling = ['annular MLR F475W','annular MLR F814W','annular MLR F160W']
    blabeling = ['single MLR F475W','single MLR F814W','single MLR F160W']
    
# up to this point, jr_comp and jr_overlay were identical

    # we change the radii into kpc (do we need a new thing for this, or can i use the old?)

    # plot specific anular mass (from aMLR), using only data from 814 really.
    ax = fig.add_subplot(2,1,1)
    ax.plot(kpc_radius, aMass[1], acolors[1], marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, aMass[1], acolors[1])
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.title(galaxies[0] + 'Mass vs. Radius, annular M/L ratios',fontsize=16)
    plt.tight_layout()
    legend = ax.legend(loc='upper right')
    
    # plot broad annular mass (light mass) from bMLR
    bx = fig.add_subplot(2,1,1)
    bx.plot(kpc_radius, bMass[1], 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bMass[1], 'yellowgreen')
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.tight_layout()
    plt.title(galaxies[0] + ' Mass vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')
    
# plot mass density per area vs radius
    # calculate mass surface density (MSD) in solarmasses /kpc^2 for aMass, bMass
    aMSD = aMass[1]/areaKpc
    bMSD = bMass[1]/areaKpc
    
    # now plotting aMass surface density vs radius [solarmasses / kpc^2]
    ax = fig.add_subplot(2,1,2)
    ax.plot(kpc_radius, aMSD, 'g--', marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, aMSD, 'g--')
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Surface Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title(galaxies[0]+ 'Mass Surface Density vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')
    
    # same for bMass
    bx = fig.add_subplot(2,1,2)
    bx.plot(kpc_radius, bMSD, 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bMSD, 'yellowgreen')
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Surface Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title(galaxies[0] + ' Mass Surface Density vs. Radius, single and annular MLRs',fontsize=15)
    legend = bx.legend(loc='upper right')
    
    # calculate and print total aMass, bMass for annular and single MLR calculations
    # using colorUV, aColorUV
    print('total mass, annular', np.sum(aMass[1])
    print('total mass, broad', np.sum(bMass[1])

    pdf.savefig()
    plt.close()

    os.system('open %s &' % 'sg_comp1_' + galaxies[0] + '.pdf')



