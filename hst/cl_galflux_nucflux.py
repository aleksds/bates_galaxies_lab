#Compile measurements of total galaxy flux (based on photometry) and nuclear flux (estimated from GALFIT or GALFITM models) in F814W and F475W for each galaxy.
#Make a plot of m_475W - m_814W vs m_814W with two points for each galaxy: one corresponding to the total galaxy flux and one corresponding to the nuclear component.
import os
import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# define the directory that contains the images and lists of needed stuff
dir = os.environ['HSTDIR']
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
wavelengths = ['F475W','F814W']
xcen = [3629, 3934, 3387, 3477, 3573, 3803, 3884, 4147, 3787, 4175, 3567, 4067]
ycen = [4154, 4137, 3503, 3405, 3339, 4170, 4165, 3922, 4186, 3827, 3436, 4054]
width, pixr = 81, 40
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
galenergy = np.zeros([len(galaxies),len(wavelengths), width, width])

fnudfsf = np.zeros([12])
expdfsf = np.zeros([12])

fnudef = np.zeros([12])
expdef = np.zeros([12])

fnumfsf = np.zeros([12])
expmfsf = np.zeros([12])

fnumef = np.zeros([12])
expmef = np.zeros([12])

totdataenergyfsf = np.zeros([12])
totdataenergyef = np.zeros([12])
totmodelenergyfsf = np.zeros([12])
totmodelenergyef = np.zeros([12])

modelenergyfsf = [[],[],[],[],[],[],[],[],[],[],[],[]]
modelenergyef = [[],[],[],[],[],[],[],[],[],[],[],[]]

foldername = [['475coarsefinepsf'],['814coarsefinepsf']]
fitsfilename = [['/J0826_F475W_coarse.fits', '/J0826_F814W_coarse.fits'], ['/J0901_F475W_coarse.fits', '/J0901_F814W_coarse.fits'], ['/J0905_F475W_coarse.fits', '/J0905_F814W_coarse.fits'], ['/J0944_F475W_coarse.fits', '/J0944_F814W_coarse.fits'], ['/J1107_F475W_coarse.fits', '/J1107_F814W_coarse.fits'], ['/J1219_F475W_coarse.fits', '/J1219_F814W_coarse.fits'], ['/J1341_F475W_coarse.fits', '/J1341_F814W_coarse.fits'], ['/J1506_F475W_coarse.fits', '/J1506_F814W_coarse.fits'], ['/J1558_F475W_coarse.fits', '/J1558_F814W_coarse.fits'], ['/J1613_F475W_coarse.fits', '/J1613_F814W_coarse.fits'], ['/J2116_F475W_coarse.fits', '/J2116_F814W_coarse.fits'], ['/J2140_F475W_coarse.fits', '/J2140_F814W_coarse.fits']]

#read in flux from data and energy from galfit nuclear model
for w in range (0, 12):
    for i in range (0, 2):
            
        file = glob.glob(dir+galaxies[w]+'*/coarse/'+wavelengths[i]+'/final*sci.fits')
        hdu = fits.open(file[0]) 
        data[i], header[i] = hdu[0].data, hdu[0].header
        if i == 0:
            fnudfsf[w] = header[i]['PHOTFNU']
            expdfsf[w] = header[i]['EXPTIME']
        if i == 1:
            fnudef[w] = header[i]['PHOTFNU']
            expdef[w] = header[i]['EXPTIME']
        for j in range(0,width):
            for k in range(0,width):
                galenergy[w,i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1]))

        galfile = glob.glob('/Volumes/physics/linux-lab/data/galfit/'+foldername[i][0]+fitsfilename[w][i])
        mc = fits.open(galfile[0])
        mc_model, mc_header = mc[2].data, mc[3].header
        
        if i == 0:
            modelenergyfsf[w] = mc_model
            fnumfsf[w] = header[i]['PHOTFNU']
            expmfsf[w] = header[i]['EXPTIME']
        if i == 1:
            modelenergyef[w] = mc_model
            fnumef[w] = header[i]['PHOTFNU']
            expmef[w] = header[i]['EXPTIME']
        
#carving out the 81 by 81 square centered on the centroid for the model energyes
shavenmodelenergyfsf = np.zeros([12,width,width])
shavenmodelenergyef = np.zeros([12,width,width])
for w in range (0,12):
    for p in range(width):
        for q in range(width):
            shavenmodelenergyfsf[w][p][q] = modelenergyfsf[w][161+p][161+q]
            shavenmodelenergyef[w][p][q] = modelenergyef[w][161+p][161+q]
            
#summing up energyes to get total energy for each galaxy in the two filters from the model
for w in range(0,12):
    totmodelenergyfsf[w] = np.sum(shavenmodelenergyfsf[w][:][:])
    totmodelenergyef[w] = np.sum(shavenmodelenergyef[w][:][:])
        
#summing up energyes to get total energy for each galaxy in the two filters from the data
for w in range(0,12):
    for i in range(0,2):
        if i == 0:
            totdataenergyfsf[w] = np.sum(galenergy[w][i][:][:])
#            totmodelenergyfsf[w] = np.sum(shavenmodelenergyfsf[w][i])
        if i == 1:
            totdataenergyef[w] = np.sum(galenergy[w][i][:][:])

#above shit aint flux its energy! Now dividing by exposure time to get power. Then will multiply by photfnu to get total flux in janskys which will allow me to calculate magnitudes! YEYEYE
fluxdatafsf = (totdataenergyfsf/expdfsf)*fnudfsf
fluxdataef = (totdataenergyef/expdef)*fnudef
fluxmodelfsf = (totmodelenergyfsf/expmfsf)*fnumfsf
fluxmodelef = (totmodelenergyef/expmef)*fnumef

#MAGNITUDES
magdatafsf = -2.5*np.log10(fluxdatafsf)
magdataef = -2.5*np.log10(fluxdataef)
magmodelfsf = -2.5*np.log10(fluxmodelfsf)
magmodelef = -2.5*np.log10(fluxmodelef)

with PdfPages('Color Plot.pdf') as pdf:   
    plt.figure()
    plt.scatter(magdataef,magdatafsf-magdataef, label='Data')
    plt.scatter(magmodelef,magmodelfsf-magmodelef, label='Model')
  #  plt.xlim(17,19)
   # plt.ylim(0,22)
    plt.xlabel('Mag 814')
    plt.ylabel('Mag 475 - Mag 814')
    plt.title('Color Plot')
    plt.legend(loc='upper right')
    pdf.savefig()
    plt.close()
#leaving off: mc_model for 475 is WHACK
