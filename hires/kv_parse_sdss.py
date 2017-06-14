#write code that will look into sdss files and write loop for each galaxy figure out each ra and dec
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
dir = os.environ['SDSSDIR']
file = glob.glob(dir+'*.fits')
#read in the relevant sdss files
for i in range(0, len(file)):
    hdulist = fits.open(file[i])
    plug_ra = hdulist[0].header['PLUG_RA']
    plug_dec = hdulist[0].header['PLUG_DEC']
    print(plug_ra) #prints ra values for each of galaxies
    print(plug_dec)#prints dec values for each galaxy 
    c = SkyCoord(ra=plug_ra*u.degree, dec=plug_dec*u.degree)
    test = c.to_string('hmsdms')#this variable gives the sky location of each galaxy
    print(file[i], test)
