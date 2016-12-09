# Joshua Rines
# 20161128
#
#the goal of this script is to take our flux values, whether it be in annulus form or total flux form, and calculate a few things:
# 1. luminosity distance using Ned Wright's calculator (use calc. online, then write code for it later)
# 2. measure total flux in galaxy
# 3. compute luminosity using Hogg equation 24
# 4. compute mass-to-light ratio using Bell & de Jong Table 3
    #using 7 a,b value sets to get a best average for masses
# 5. compute mass of galaxy and annular-based mass
    #Msic and Msrc calculations based on spatially-integrated and spatially-reolved MLR
# 6. for galaxy J0826, measure flux, color, luminosity, mass-to-light ratio and mass for each annulus
# 7. make plots of spatially-resolved mass vs. radius and spatially-integrated mass vs. radius for J0826
    #include dividing out by the area of the annnulus, and making conversion from pixels to kpc
    #also only considering the 'green' filter (F814W) here so we can make better use of the best value MLRs
    #overlay of the Msrc and Msic plots for mass vs. radius and mass/area vs. radius (kpc)

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

# specify the position of the science target and the size of the region around the science target to consider
xcen = 3628.7
ycen = 4153.8
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
with PdfPages('jr_overlays_J0826.pdf') as pdf:

    fig = plt.figure()

    collection = ['F475W','F814W','F160W']

    flux = np.zeros([len(collection),len(radii)])
    subflux = np.zeros([len(collection),len(radii)])

    for i in range (0, len(collection)):
        
        # read in the images
        file = glob.glob(dir+'J0826_final_'+collection[i]+'*sci.fits')
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

    #luminosity distance (in cm) for J0826, z = 0.603
    LdJ0826 = 3551.4*3.08568*10**24

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
    LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0826**2)
    LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0826**2)
    LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0826**2)

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

    mass = (M_BV_Bk, M_BV_Vk, M_BV_Jk)

    #calculating best values and uncertainties
    mass = (M_BV_Bk, M_BV_Vk, M_BV_Jk)
    Msic_475_BV = np.mean(mass[0])
    Msic_814_BV = np.mean(mass[1])
    Msic_160_BV = np.mean(mass[2])
    
    Msic_475_BV_std = np.std(mass[0])
    Msic_814_BV_std = np.std(mass[1])
    Msic_160_BV_std = np.std(mass[2])
    
    print('Msic,475W,B-V', Msic_475_BV/1e11)
    print('Msic,475W,B-V std', Msic_475_BV_std/1e11)
    
    print('Msic,814W,B-V', Msic_814_BV/1e11)
    print('Msic,814W,B-V std', Msic_814_BV_std/1e11)
    
    print('Msic,160W,B-V', Msic_160_BV/1e11)
    print('Msic,160W,B-V std', Msic_160_BV_std/1e11)

    #determining M/L ratio using Table 1 of Bell & de Jong
    #MLR_BV_B = 10**(-.994+(1.804*colorUV))
    #MLR_BV_V = 10**(-.734+(1.404*colorUV))
    #MLR_BV_J = 10**(-0.621+(0.794*colorUV))
    #MLR_VJ_B = 10**(-1.903+(1.138*colorVJ))
    #MLR_VJ_V = 10**(-1.477+0.905*(colorVJ))
    #MLR_VJ_J = 10**(-1.029+(.505*colorVJ))

    #calculating nu_e * L_nu_e luminosity in erg/s units from Hogg eq (24)
    #c = 299792458
    #LnuNu475 = (c/(475*10**-9))*tflux475*(4*math.pi*LdJ0826**2)
    #LnuNu814 = (c/(814*10**-9))*tflux814*(4*math.pi*LdJ0826**2)
    #LnuNu160 = (c/(1600*10**-9))*tflux160*(4*math.pi*LdJ0826**2)

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

    #after discussing with Aleks, we determined that we do want an individualized MLR for each annulus, so I have done that by replacing the 'acolorUV' below. This gives an array with size 40 with which we can calculate the aLunuNu, and later on the mass in each ring.
    aMLR_BV_B = 10**(-.994+(1.804*acolorUV))
    aMLR_BV_V = 10**(-.734+(1.404*acolorUV))
    aMLR_BV_J = 10**(-.621+(.794*acolorUV))

    #setting up to calculate Msrc_814_BV_ab0-6
    for k in range(len(Va)):
        #aMLR_BV_Vk_ab0 = np.zeros(40)
        aMLR_BV_Vk_ab0 = 10**(Va[0]+(Vb[0]*acolorUV))
        aMLR_BV_Vk_ab1 = 10**(Va[1]+(Vb[1]*acolorUV))
        aMLR_BV_Vk_ab2 = 10**(Va[2]+(Vb[2]*acolorUV))
        aMLR_BV_Vk_ab3 = 10**(Va[3]+(Vb[3]*acolorUV))
        aMLR_BV_Vk_ab4 = 10**(Va[4]+(Vb[4]*acolorUV))
        aMLR_BV_Vk_ab5 = 10**(Va[5]+(Vb[5]*acolorUV))
        aMLR_BV_Vk_ab6 = 10**(Va[6]+(Vb[6]*acolorUV))
        aMLR_BV_Vk_ab = (aMLR_BV_Vk_ab0,aMLR_BV_Vk_ab1,aMLR_BV_Vk_ab2,aMLR_BV_Vk_ab3,aMLR_BV_Vk_ab4,aMLR_BV_Vk_ab5,aMLR_BV_Vk_ab6)

    #calculation of annular MLR based on one MLR for the entire galaxy, since we plan to overlay polots of both. this will be denoted 'bMLR...'
    bMLR_BV_B = 10**(-.994+(1.804*colorUV))
    bMLR_BV_V = 10**(-.734+(1.404*colorUV))
    bMLR_BV_J = 10**(-.621+(.794*colorUV))

    #calculating nu_e * L_nu_e luminosity in erg/s units for each annulus from Hogg eq (24)
    c = 299792458
    aLnuNu475 = (c/(475*10**-9))*aflux475*(4*math.pi*LdJ0826**2)
    aLnuNu814 = (c/(814*10**-9))*aflux814*(4*math.pi*LdJ0826**2)
    aLnuNu160 = (c/(1600*10**-9))*aflux160*(4*math.pi*LdJ0826**2)

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

    #calculating the Msrc values for each annulus
    aMsrc_814_BV_ab0 = aLsol814*aMLR_BV_Vk_ab[0]
    aMsrc_814_BV_ab1 = aLsol814*aMLR_BV_Vk_ab[1]
    aMsrc_814_BV_ab2 = aLsol814*aMLR_BV_Vk_ab[2]
    aMsrc_814_BV_ab3 = aLsol814*aMLR_BV_Vk_ab[3]
    aMsrc_814_BV_ab4 = aLsol814*aMLR_BV_Vk_ab[4]
    aMsrc_814_BV_ab5 = aLsol814*aMLR_BV_Vk_ab[5]
    aMsrc_814_BV_ab6 = aLsol814*aMLR_BV_Vk_ab[6]
    aMsrc_814_BV_ab = (aMsrc_814_BV_ab0,aMsrc_814_BV_ab1,aMsrc_814_BV_ab2,aMsrc_814_BV_ab3,aMsrc_814_BV_ab4,aMsrc_814_BV_ab5,aMsrc_814_BV_ab6)

    #calculating the 7 Msrc values, one from each set of ab values
    Msrc_814_BV_ab0 = np.sum(aMsrc_814_BV_ab[0])
    Msrc_814_BV_ab1 = np.sum(aMsrc_814_BV_ab[1])
    Msrc_814_BV_ab2 = np.sum(aMsrc_814_BV_ab[2])
    Msrc_814_BV_ab3 = np.sum(aMsrc_814_BV_ab[3])
    Msrc_814_BV_ab4 = np.sum(aMsrc_814_BV_ab[4])
    Msrc_814_BV_ab5 = np.sum(aMsrc_814_BV_ab[5])
    Msrc_814_BV_ab6 = np.sum(aMsrc_814_BV_ab[6])
    Msrc_814_BV_ab = (Msrc_814_BV_ab0,Msrc_814_BV_ab1,Msrc_814_BV_ab2,Msrc_814_BV_ab3,Msrc_814_BV_ab4,Msrc_814_BV_ab5,Msrc_814_BV_ab6)

    #best value for each annulus
    bestval_annular_Msrc = np.zeros(40)
    for j in range(len(radii)):
        bestval_annular_Msrc[j] = (aMsrc_814_BV_ab0[j]+aMsrc_814_BV_ab1[j]+aMsrc_814_BV_ab2[j]+aMsrc_814_BV_ab3[j]+aMsrc_814_BV_ab4[j]+aMsrc_814_BV_ab5[j]+aMsrc_814_BV_ab6[j])/7

    #best value and std, printed
    Msrc_814_BV = np.mean(Msrc_814_BV_ab)
    Msrc_814_BV_std = np.std(Msrc_814_BV_ab)
    print('Msrc,814W,B-V',Msrc_814_BV/1e11)
    print('Msrc,814W,B-V std',Msrc_814_BV_std/1e11)

    #calculate mass associated with each annulus in solar units, based on one MLR estimate for the entire galaxy
    bM_BV_B = aLsol475*bMLR_BV_B
    bM_BV_V = aLsol814*bMLR_BV_V
    bM_BV_J = aLsol160*bMLR_BV_J
    bmass = (bM_BV_B,bM_BV_V,bM_BV_J)

    #best values
    annular_Msic_814_BV_ab0 = aLsol814*MLR_BV_Vk[0]
    annular_Msic_814_BV_ab1 = aLsol814*MLR_BV_Vk[1]
    annular_Msic_814_BV_ab2 = aLsol814*MLR_BV_Vk[2]
    annular_Msic_814_BV_ab3 = aLsol814*MLR_BV_Vk[3]
    annular_Msic_814_BV_ab4 = aLsol814*MLR_BV_Vk[4]
    annular_Msic_814_BV_ab5 = aLsol814*MLR_BV_Vk[5]
    annular_Msic_814_BV_ab6 = aLsol814*MLR_BV_Vk[6]
    annular_Msic_814_BV_ab = (annular_Msic_814_BV_ab0,annular_Msic_814_BV_ab1,annular_Msic_814_BV_ab2,annular_Msic_814_BV_ab3,annular_Msic_814_BV_ab4,annular_Msic_814_BV_ab5,annular_Msic_814_BV_ab6)
    
    #best value for each annulus
    bestval_annular_Msic = np.zeros(40)
    for j in range(len(radii)):
        bestval_annular_Msic = (annular_Msic_814_BV_ab[0]+annular_Msic_814_BV_ab[1]+annular_Msic_814_BV_ab[2]+annular_Msic_814_BV_ab[3]+annular_Msic_814_BV_ab[4]+annular_Msic_814_BV_ab[5]+annular_Msic_814_BV_ab[6])/7

    #getting seven values for total mass of the galaxy
    Msic_814_BV_ab0 = np.sum(annular_Msic_814_BV_ab0)
    Msic_814_BV_ab1 = np.sum(annular_Msic_814_BV_ab1)
    Msic_814_BV_ab2 = np.sum(annular_Msic_814_BV_ab2)
    Msic_814_BV_ab3 = np.sum(annular_Msic_814_BV_ab3)
    Msic_814_BV_ab4 = np.sum(annular_Msic_814_BV_ab4)
    Msic_814_BV_ab5 = np.sum(annular_Msic_814_BV_ab5)
    Msic_814_BV_ab6 = np.sum(annular_Msic_814_BV_ab6)
    Msic_814_BV_ab = (Msic_814_BV_ab0,Msic_814_BV_ab1,Msic_814_BV_ab2,Msic_814_BV_ab3,Msic_814_BV_ab4,Msic_814_BV_ab5,Msic_814_BV_ab6)
    #NOTE: taking the mean of Msic_814_BV_ab will give you the same thing as Msic_814_BV

    #making various calculations, including std, best value, etc.
    #using a different method, but I'll just comment this out in case we need it later
    #ave_mass475 = np.mean(mass[0])
    #ave_mass814 = np.mean(mass[1])
    #ave_mass160 = np.mean(mass[2])
    #ave_mass = (ave_mass475,ave_mass814,ave_mass160)
    
    #std_mass475 = np.std(mass[0])
    #std_mass814 = np.std(mass[1])
    #std_mass160 = np.std(mass[2])
    #std_mass = (std_mass475,std_mass814,std_mass160)

    #step1 = ave_mass[0]/np.square(std_mass[0]) + ave_mass[1]/np.square(std_mass[1]) + ave_mass[2]/np.square(std_mass[2])
    #step2 = 1/np.square(std_mass[0]) + 1/np.square(std_mass[1]) + 1/np.square(std_mass[2])

    #best_value = step1/step2
    #print('J0826 total mass best value:', best_value/1e11)

    #uncert = (1/np.square(std_mass[0]) + 1/np.square(std_mass[1]) + 1/np.square(std_mass[2]))**-.5
    #print('J0826 total mass uncertainty:', uncert/1e11)

    #for an individual annulus
    #for j in range(0,len(collection)+1):
        #amass_0 = amass[j][0]+amass[j+1][0]+amass[j+2][0]

    #plotting mass vs radius
    acolors = ['b--','g--','r--']
    bcolors = ['b', 'g', 'r']
    adot = ['b','g','r']
    bdot = ['bo','go','ro']
    alabeling = ['annular MLR F475W','M_SRC F814W','annular MLR F160W']
    blabeling = ['single MLR F475W','M_SIC F814W','single MLR F160W']

    #putting radius into kpc
    kpc_radius = radii*(0.05)*(6.701)

    #plotting the specific annular (specific aMLR) mass
    #for k in range(0,len(acolors)):
    ax = fig.add_subplot(2,1,1)
    ax.plot(kpc_radius, bestval_annular_Msrc, acolors[1], marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, bestval_annular_Msrc, acolors[1])
    #plt.plot(np.unique(radii), np.poly1d(np.polyfit(radii, amass[k], 192))(np.unique(radii)),bcolors[k], label=str(alabeling[k]))
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.title('J0826 Mass vs. Radius, annular M/L ratios',fontsize=16)
    plt.tight_layout()
    legend = ax.legend(loc='upper right')
        
    #plotting the broad (single bMLR) annular mass, the 'light mass' as I call it
    #for k in range(0,len(bcolors)):
    bx = fig.add_subplot(2,1,1)
    bx.plot(kpc_radius, bestval_annular_Msic, 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bestval_annular_Msic, 'yellowgreen')
    #plt.plot(np.unique(radii), np.poly1d(np.polyfit(radii, bmass[k], 192))(np.unique(radii)),bcolors[k], label=str(alabeling[k]))
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass (solar masses)',fontsize=14)
    plt.tight_layout()
    plt.title('J0826 Mass vs. Radius, M_SRC and M_SIC',fontsize=15)
    legend = bx.legend(loc='upper right')

    
    #plotting the mass/area vs radius, converting radius units from pixels to kpc to get mass surface density units
    
    #first creating an array with areas of shells in proper units of kpc
    kpc_area = np.zeros(len(area))
    bestval_annular_Msrc_ovr_area = np.zeros(len(area))
    bestval_annular_Msic_ovr_area = np.zeros(len(area))
    for j in range(0,len(area)):
       kpc_area[j] = area[j]*(0.05**2)*(6.701**2)
       
    #now calculating mass/area (mass surface density) in units of solar masses/kpc for bestval_annular_Msrc and bestval_annular_Msic and also radius in kpc units
    bestval_annular_Msrc_ovr_area = amass[1]/kpc_area
    bestval_annular_Msic_ovr_area = bmass[1]/kpc_area
    kpc_radius = radii*(0.05)*(6.701)
    
    #now plotting bestval_annular_Msrc_ovr_area and bestval_annular_Msic_ovr_area vs radius in kpc
    ax = fig.add_subplot(2,1,2)
    ax.plot(kpc_radius, bestval_annular_Msrc_ovr_area, 'g--', marker='s', label=str(alabeling[1]))
    ax.plot(kpc_radius, bestval_annular_Msrc_ovr_area, 'g--')
    #plt.plot(np.unique(radii), np.poly1d(np.polyfit(radii, bmass[k], 192))(np.unique(radii)),bcolors[k], label=str(alabeling[k]))
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title('J0826 Mass Density vs. Radius, M_SRC and M_SIC',fontsize=15)
    legend = bx.legend(loc='upper right')

    #plotting bestval_annular_Msic_ovr_area vs radius in kpc
    bx = fig.add_subplot(2,1,2)
    bx.plot(kpc_radius, bestval_annular_Msic_ovr_area, 'yellowgreen', marker='o', label=str(blabeling[1]))
    bx.plot(kpc_radius, bestval_annular_Msic_ovr_area, 'yellowgreen')
    #plt.plot(np.unique(radii), np.poly1d(np.polyfit(radii, bmass[k], 192))(np.unique(radii)),bcolors[k], label=str(alabeling[k]))
    plt.xlabel('Radius (kpc)',fontsize=14)
    plt.ylabel('Mass Density (M_sol/area)',fontsize=14)
    plt.tight_layout()
    plt.title('J0826 Mass Density vs. Radius, M_SRC and M_SIC',fontsize=15)
    legend = bx.legend(loc='upper right')
    
    
    #calculating the area under the curve for amass (F814W filter only)
    #ended up not needing to do this, but leaving it in the code
    #atrap = np.zeros([len(radii)-1])
    #for j in range(0,39):
        #atrap[j] = .5*(amass[1][j]+amass[1][j+1])
    #aintegral = np.sum(atrap)
    #half_aintegral = aintegral / 2

    #calculating the area under the curve for bmass (F814 filter only)
    #ended up not needing to do this, but leaving it in the code
    #btrap = np.zeros([len(radii)-1])
    #for j in range(0,39):
        #btrap[j] = .5*(bmass[1][j]+bmass[1][j+1])
    #bintegral = np.sum(btrap)
    #half_bintegral = bintegral / 2

    #calculating total mass (amass) for annular MLR (814 filter only)
    total_annular_Msrc_F814W = np.sum(bestval_annular_Msrc)
    print('Msrc,814,BV total', total_annular_Msrc_F814W/1e11)

    #calculating total mass (bmass) for single MLR (814 filter only)
    total_singular_Msic_F814W = np.sum(bestval_annular_Msic)
    print('Msic,814,BV total', total_singular_Msic_F814W/1e11)

    #calculating %amass and %bmass in first 5 annuli
    Msrc_first_5 = np.sum(bestval_annular_Msrc[0:4])
    pct_Msrc_first_5 = Msrc_first_5/total_annular_Msrc_F814W*100
    Msic_first_5 = np.sum(bestval_annular_Msrc[0:4])
    pct_Msic_first_5 = Msic_first_5/total_singular_Msic_F814W*100
    print('% Msrc first 5', pct_Msrc_first_5)
    print('% Msic first 5', pct_Msic_first_5)


    pdf.savefig()
    plt.close()

    
    os.system('open %s &' % 'jr_overlays_J0826.pdf')
