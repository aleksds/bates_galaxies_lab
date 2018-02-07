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

#PLOT CODE STELLAR MASS VS R_E

#masses = stellar masses in the order: J0826, J0901, J0905, J0944, J1107, J1219, J1341, J1506, J1558, J1613, J2116, J2140
masses = [11.3883, 11.5977, 11.2696, 11.1481, 11.0063, 11.2015, 10.8927, 10.9912, 11.1112, 11.4324, 11.2846, 11.1891]

filename = 'cl_masskpcnewpsfplot.pdf'
with PdfPages(filename) as pdf:
    fig = plt.figure()
    plt.title('Stellar Mass vs. Effective Radius' ,fontsize=20)
    plt.xlabel('Stellar Mass (M☉)', fontsize=17)
    plt.ylabel('Effective Radius (kpc)', fontsize=17)
    plt.yscale('log')
    
    plt.xlim((10.8,11.7))
    plt.xticks([10.9,11.1,11.3,11.5,11.7])
    plt.ylim((0.05,2))
    #plt.yticks([0.1,0.2,1])
    plt.tick_params(axis='x',labelsize=15)
    plt.tick_params(axis='y',labelsize=15)

    colors = ['purple','blue','red','orange']
    markers = ['x','D','o','>']
    
    for c in range(12):
        for i in range(4):
            plt.scatter(masses[c],storage[c][i], s=30, marker = markers[i], c=colors[i])
    a = plt.scatter(masses[1],storage[1][0], s=30, marker = markers[0], c=colors[0])
    b = plt.scatter(masses[1],storage[1][1], s=30, marker = markers[1], c=colors[1])
    c = plt.scatter(masses[1],storage[1][2], s=30, marker = markers[2], c=colors[2])
    d = plt.scatter(masses[1],storage[1][3], s=30, marker = markers[3], c=colors[3])
    plt.legend((a,b,c,d),
           ('475 coarse with fine psf', '475 fine', '814 coarse with fine psf', '814 fine'),
           scatterpoints=1,
           loc='upper left',
           ncol=2,
           fontsize=8)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
os.system('open %s &' % filename)


#475 coarse with fine psf, 475 fine, 814 coarse with fine psf, 814 fine
