#Cristopher Thopmson
#Date: 2019/04/25
#The purpose of this code is to pull and print the chi values for the simulatneous coarse 
#values produced by the input_creator_V5.py. This code takes the values for the seven 
#different band files for each of the 12 galaxies.

import sys
import os
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from photutils import CircularAperture
from photutils import aperture_photometry

dir = sys.argv[1]
files = glob.glob(dir+'/*.band')
images = glob.glob(dir+'/*output.fits')
chi = np.zeros(len(files))
file = "All_the_Chi_Values.dat"
file_object = open('All_the_Chi_Values.dat','a+')
file_object.write("Chi Values From Galfit\n")
chi_from_our_calculations = np.zeros(len(files))
chi_from_our_calculations_with_dof = np.zeros(len(files))
chi = np.zeros(len(files))
chi2 = np.zeros(len(files))
for i in range (0, len(files)):
    print(images[i])
    hdu = fits.open(images[i])
    # F814W corresponds roughly to rest-frame V-band
    vdat, vdat_head = hdu[0].data, hdu[0].header
    vunc, vunc_head = hdu[12].data, hdu[12].header
    vmod, vmod_head = hdu[3].data, hdu[3].header
    vres, vres_head = hdu[6].data, hdu[6].header

    # F475W corresponds roughly to rest-frame U-band
    udat, udat_head = hdu[1].data, hdu[1].header
    uunc, uunc_head = hdu[13].data, hdu[13].header
    umod, umod_head = hdu[4].data, hdu[4].header
    ures, ures_head = hdu[7].data, hdu[7].header

    # F160W corresponds roughly to rest-frame J-band
    jdat, jdat_head = hdu[2].data, hdu[2].header
    junc, junc_head = hdu[14].data, hdu[14].header
    jmod, jmod_head = hdu[5].data, hdu[5].header
    jres, jres_head = hdu[8].data, hdu[8].header
    chi_from_our_calculations[i] = np.sum((vres/vunc)**2) + np.sum((ures/uunc)**2) + np.sum((jres/junc)**2)
    chi_from_our_calculations_with_dof[i] = np.sum((vres/vunc)**2)/((101**2)*3)
    print(files[i])
    with open(files[i]) as f:
        content = f.readlines()
        print(content[3][15:19])
        chi2[i] = np.float(content[3][15:19])
        print(content[3][31:40])
        chi[i] = np.float(content[3][31:40])
    print(chi_from_our_calculations[i])
    print(chi_from_our_calculations_with_dof[i])
    file_object.write(files[i]+' '+str(chi[i])+" "+str(chi2[i])+" "+str(chi_from_our_calculations[i])+'\n')
file_object.close
# Up to this part, it creates a code which works and produces all of the chi values in a colum
# and a name is printed out for each of the data that are given by the galfit file. 
# Upward and onward to coming up with code for the calculated values of chi

