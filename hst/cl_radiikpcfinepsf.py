import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from xlrd import open_workbook
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

#redshifts in order:J0826,J0901,J0905,J0944,J1107,J1219,J1341,J1506,J1558,J1613,J2116,J2140
redshifts = [0.603,0.459,0.712,0.514,0.467,0.451,0.451,0.658,0.608,0.402,0.449,0.728,0.752]

#making mad lists (content from GALFIT)
#galaxy order: J0826,J0901,J0905,J0944,J1107,J1219,J1341,J1506,J1558,J1613,J2116,J2140
#radii in order:475 coarse with fine psf, 475 fine, 814 coarse with fine psf, 814 fine
#some are zeros because i need to run galfit on certain images and some are 0.01 because i need to fix initial parameters and run galfit again
galaxyradiipix = [[0.25,0.46,0.4,0.87],[0.6,0.81,0.51,1.06],[0.24,0.41,0.23,0.45],[0.51,0.47,0.24,0.6],[0.21,0.49,0.61,1.26],[0.57,.9,0.83,1.66],[0.59,0.99,0.28,0.56],[0.22,0.4,0.31,0.63],[0.64,1.22,3.07,5.73],[2.08,3.91,3.29,6.24],[0.34,0.63,0.62,1.15],[0.25,0.22,0.52,0.94]]

J0826radiikpc = np.zeros(4)
J0901radiikpc = np.zeros(4)
J0905radiikpc = np.zeros(4)
J0944radiikpc = np.zeros(4)
J1107radiikpc = np.zeros(4)
J1219radiikpc = np.zeros(4)
J1341radiikpc = np.zeros(4)
J1506radiikpc = np.zeros(4)
J1558radiikpc = np.zeros(4)
J1613radiikpc = np.zeros(4)
J2116radiikpc = np.zeros(4)
J2140radiikpc = np.zeros(4)
storage = [J0826radiikpc,J0901radiikpc,J0905radiikpc,J0944radiikpc,J1107radiikpc,J1219radiikpc,J1341radiikpc,J1506radiikpc,J1558radiikpc,J1613radiikpc,J2116radiikpc,J2140radiikpc]

#functions to convert from pix to arcseconds
def pixtoarcseccoarse(pix):
    arcsecradius = pix*(0.05)
    return arcsecradius
def pixtoarcsecfine(pix):
    arcsecradius = pix*(0.025)
    return arcsecradius

#loop to convert the radius from pix to kpc for a single galaxy    
for w in range(0,12):
    q = galaxyradiipix[w]
    arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshifts[w])
    for i in range(0, 4):
        if i == 0 or i == 2:
            storage[w][i] = round(pixtoarcseccoarse(q[i])/arcsecperkpc.value, 3)
        if i == 1 or i == 3:
            storage[w][i] = round(pixtoarcsecfine(q[i])/arcsecperkpc.value, 3)




#extracting the 2 re that we care about (475 coarse with fine psf, 814 coarse with fine psf)
J0826kpc = np.zeros(2)
J0901kpc = np.zeros(2)
J0905kpc = np.zeros(2)
J0944kpc = np.zeros(2)
J1107kpc = np.zeros(2)
J1219kpc = np.zeros(2)
J1341kpc = np.zeros(2)
J1506kpc = np.zeros(2)
J1558kpc = np.zeros(2)
J1613kpc = np.zeros(2)
J2116kpc = np.zeros(2)
J2140kpc = np.zeros(2)
for w in range(0,4):
    if w == 0:
        J0826kpc[w] = J0826radiikpc[w]
        J0901kpc[w] = J0901radiikpc[w]
        J0905kpc[w] = J0905radiikpc[w]
        J0944kpc[w] = J0944radiikpc[w]
        J1107kpc[w] = J1107radiikpc[w]
        J1219kpc[w] = J1219radiikpc[w]
        J1341kpc[w] = J1341radiikpc[w]
        J1506kpc[w] = J1506radiikpc[w]
        J1558kpc[w] = J1558radiikpc[w]
        J1613kpc[w] = J1613radiikpc[w]
        J2116kpc[w] = J2116radiikpc[w]
        J2140kpc[w] = J2140radiikpc[w]
    if w == 2:
        J0826kpc[w-1] = J0826radiikpc[w]
        J0901kpc[w-1] = J0901radiikpc[w]
        J0905kpc[w-1] = J0905radiikpc[w]
        J0944kpc[w-1] = J0944radiikpc[w]
        J1107kpc[w-1] = J1107radiikpc[w]
        J1219kpc[w-1] = J1219radiikpc[w]
        J1341kpc[w-1] = J1341radiikpc[w]
        J1506kpc[w-1] = J1506radiikpc[w]
        J1558kpc[w-1] = J1558radiikpc[w]
        J1613kpc[w-1] = J1613radiikpc[w]
        J2116kpc[w-1] = J2116radiikpc[w]
        J2140kpc[w-1] = J2140radiikpc[w]
  
