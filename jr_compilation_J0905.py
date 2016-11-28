# Joshua Rines
# 20161122
#
#the goal of this script is to take our flux values, whether it be in annulus form or total flux form, and calculate a few things:
# 1. luminosity distance using Ned Wright's calculator (use calc. online, then write code for it later)
# 2. measure total flux in galaxy
# 3. compute luminosity using Hogg equation 24
# 4. compute mass-to-light ratio using Bell & de Jong Table 1
# 5. compute mass of galaxy
# 6. for galaxy J0905, measure flux, color, luminosity, mass-to-light ratio and mass for each annulus
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
with PdfPages('jr_compilation_J0905.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir+'J0905_final_'+collection[i]+'*sci.fits')
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
        

    #calculating galaxy-wide

    #set up two-dimensional arrays for the a and b coefficients based on the luminosity and color
    #this will be in the same format ish as the table in Josh's blue notebook
    Ba = [-1.019,-1.113,-1.026,-.990,-1.110,-.994,-.888]
    Bb = [1.937,2.065,1.954,1.883,2.018,1.804,1.758]
    B_coeff = [Ba,Bb]
    Va = [-.759,-.853,-.766,-.730,-.850,-.734,-.628]
    Vb = [1.537,1.665,1.554,1.483,1.618,1.404,1.358]
    V_coeff = [Va,Vb]
    Ja = [-.540,-.658,-.527,-.514,-.659,-.621,-.550]
    Jb = [.757,.907,.741,.704,.878,.794,.801]
    J_coeff = [Ja,Jb]
    
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

    #determining M/L ratio using Table 3 of Bell & de Jong, seven coefficients for each luminosity for BV color
    MLR_BV_Bk = np.zeros([len(Ba)])
    MLR_BV_Vk = np.zeros([len(Va)])
    MLR_BV_Jk = np.zeros([len(Ja)])
    for k in range(0,len(Ba)):
        MLR_BV_Bk[k] = 10**(Ba[k]+(Bb[k]*colorUV))
        MLR_BV_Vk[k] = 10**(Va[k]+(Vb[k]*colorUV))
        MLR_BV_Jk[k] = 10**(Ja[k]+(Jb[k]*colorUV))

    #calculating nu_e * L_nu_e luminosity in erg/s units from Hogg eq (24), only three values depending on filter
    c = 299792458
    LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0905**2)
    LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0905**2)
    LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0905**2)

    #convert luminosity to solar units
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

    #determining M/L ratio using Table 1 of Bell & de Jong
    #MLR_BV_B = 10**(-.994+(1.804*colorUV))
    #MLR_BV_V = 10**(-.734+(1.404*colorUV))
    #MLR_BV_J = 10**(-0.621+(0.794*colorUV))
    #MLR_VJ_B = 10**(-1.903+(1.138*colorVJ))
    #MLR_VJ_V = 10**(-1.477+0.905*(colorVJ))
    #MLR_VJ_J = 10**(-1.029+(.505*colorVJ))

    #calculating nu_e * L_nu_e luminosity in erg/s units from Hogg eq (24)
    #c = 299792458
    #LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0905**2)
    #LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0905**2)
    #LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0905**2)

    #convert luminosity to solar units
    #Lsol475 = LnuNu475 / (3.846*10**33)
    #Lsol814 = LnuNu814 / (3.846*10**33)
    #Lsol160 = LnuNu160 / (3.846*10**33)

    #calculate mass of galaxy in solar units
    #M_BV_B = Lsol475*MLR_BV_B
    #M_BV_V = Lsol814*MLR_BV_V
    #M_BV_J = Lsol160*MLR_BV_J
    #M_VJ_B = Lsol475*MLR_VJ_B
    #M_VJ_V = Lsol814*MLR_VJ_V
    #M_VJ_J = Lsol160*MLR_VJ_J
    #print(M_BV_B/1e11, M_BV_V/1e11, M_BV_J/1e11, M_VJ_B/1e11, M_VJ_V/1e11, M_VJ_J/1e11)
    

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

    #determining M/L ratio using Table 1 of Bell & de Jong
    #need to set up a 2d array with values which will be done later, for now just using Table 1 coefficients
    #aMLR_BV_Bk = np.zeros([len(Ba)])
    #aMLR_BV_Bk = 10**(Ba[k]+(Bb[k]*colorUV))

    #I THINK I FIGURED IT OUT: I only need the ML ratio for the entire galaxy, NOT for each specific annulus...that somehow seems to mess up the total mass.  THUS: I am essentially using the same MLR as MLR_BV_X , so the 'acolor' things were not necessary I dont think...but if need be, I will put the 'acolorUV' code back in place of the 'colorUV' code below.
    
    #this is determining MLR using Table 1 coeffieicnts.  I could extend this to use all the coefficients, but for now, my plan is to calculate statistical uncertainty using the 21 values I have from finding the total mass of the galaxy and extend it to these measurements.  I could use the 3 values I have calculated below for the annuli.
    aMLR_BV_B = 10**(-.994+(1.804*colorUV))
    aMLR_BV_V = 10**(-.734+(1.404*colorUV))
    aMLR_BV_J = 10**(-.621+(.794*colorUV))

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
    aM_BV_B = aLsol475*aMLR_BV_B
    aM_BV_V = aLsol814*aMLR_BV_V
    aM_BV_J = aLsol160*aMLR_BV_J
    mass = (aM_BV_B,aM_BV_V,aM_BV_J)

    #plotting mass vs radius
    colors = ['b', 'g', 'r']
    dot = ['bo','go','ro']
    for k in range(0,len(colors)):
        ax = fig.add_subplot(1,1,1)
        ax.plot(radii,mass[k],colors[k])
        ax.plot(radii,mass[k],dot[k])
    plt.xlabel('Radius',fontsize=18)
    plt.ylabel('Mass',fontsize=18)
    plt.title('Mass vs. Radius, J0905',fontsize=20)

    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_compilation_J0905.pdf')
