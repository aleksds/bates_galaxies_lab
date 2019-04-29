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
files = glob.glob(dir+'*.band')
files = glob.glob(dir+'/*output.fits')
chi = np.zeros(len(files))
file = "All_the_Chi_Values.dat"
file_object = open('All_the_Chi_Values.dat','a+')
file_object.write("Chi Values From Galfit\n")
for i in range (0, len(files)):
    print(files[i])
    with open(files[i]) as f:
        content = f.readlines()
        print(content[3][15:19])
        chi[i] = np.float(content[3][15:19])
        vdat, vadt_head = raw_data[0].data, raw_data[0].header
        vunc, vunc_head = raw_data[12].data, raw_data[12].header
        vmod, vmod_head = raw_data[3].data, raw_data[3].header
        vres, vres_head = raw_data[6].data, raw_data[6].header
        chi_from_our_calculations[i] = np.sum((vres/vunc)**2)/((101**2*3)
        print(chi_from_our_calculations[i])
        file_object.write(str(chi[i])+ "" +str(chi_from_our_calculations[i]) '\n')
file_object.close
# Up to this part, it creates a code which works and produces all of the chi values in a colum
# and a name is printed out for each of the data that are given by the galfit file. 
# UPward and onward to coming up with code for the calculated values of chi

