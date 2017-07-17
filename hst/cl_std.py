#standard deviation of calculated centroid values

# 20170705 thoughts from Aleks
#
# (1) Could define a numpy array up front (and not append lists incrementally) with something like the following
# ngal = len(Galaxies)
# nmet = 3
# nfil = 3
# xcen = np.zeros([ngal, nmet, nfil])
# ycen = np.zeros([ngal, nmet, nfil])
# (2) then in the loop
# xcen[i] = ... output from function ...
# this would give us all the values in a way that would be useful for asking questions like:
# is the standard deviation smaller if you exclude the com values?
# is the standard deviation smaller if you exclude one of the filters?
# (3) then could calcualte standard deviation in one line of code (or whatever other quantities you wanted to look at)
# (4) then could make a plot with one or two lines of code
# (5) you could use plt.text() to put the name of the galaxy next to its point in the sigma_x vs sigma_y plot
import numpy as np
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from photutils import centroid_com, centroid_1dg, centroid_2dg

#define a function to plot "postage stamp" images
def plot_imageX(xcen,ycen):
    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    return [x1,x2,x3]

def plot_imageY(xcen,ycen):
    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    return [y1,y2,y3]
# define the directory that contains the images
dir = os.environ['HSTDIR']

# parameters for the initial guess at the galaxy centroid and the size
Galaxies = [['J0826', 3629, 4154], ['J0901', 3934, 4137], ['J0905', 3387, 3503], ['J0944', 3477., 3405.], ['J1107', 3573, 3339.], ['J1219', 3803., 4170.], ['J1341', 3884, 4165], ['J1506', 4147., 3922.], ['J1558', 3787., 4186.], ['J1613', 4175., 3827.], ['J2116', 3567, 3436], ['J2140', 4067, 4054]]

dx = 5
dy = 5

ngal = len(Galaxies)
nfil = 3
nmeth = 3
xcen = np.zeros([ngal, nmeth, nfil])
ycen = np.zeros([ngal, nmeth, nfil])

xstds = np.zeros([len(Galaxies)])
ystds = np.zeros([len(Galaxies)])

# create a PDF file for the plots
meth = ['centroid_com','centroid_1dg','centroid_2dg']
foureightone = [4,8,1]
filters = ['F475W','F814W','F160W']
colors = ['blue','green','red']

for j in range(0,len(Galaxies)):
        for i in range(0,len(filters)):
            for m in range(0,len(meth)):
                file = glob.glob(dir+Galaxies[j][0]+'_final_F'+str(foureightone[i])+'*sci.fits')
                hdu = fits.open(file[0])
                data, header = hdu[0].data, hdu[0].header
                stamp = data[round(Galaxies[j][2]-dy):round(Galaxies[j][2]+dy), round(Galaxies[j][1]-dx):round(Galaxies[j][1]+dx)]
                xcoor = plot_imageX(Galaxies[j][1],Galaxies[j][2])
                xcen[j][i][m] = xcoor[m]
                ycoor = plot_imageY(Galaxies[j][1],Galaxies[j][2])
                ycen[j][i][m] = ycoor[m]
        xstds[j] = np.std(xcen[j])
        ystds[j] = np.std(ycen[j])

with PdfPages('cl_std.pdf') as pdf:
    plt.figure()
    for q in range (0,len(Galaxies)):
        plt.scatter(xstds,ystds)
        if q==8 or q==0:
            plt.text(xstds[q]-0.0055,ystds[q]-0.006, Galaxies[q][0])
        else:
            plt.text(xstds[q]-0.0055,ystds[q]+0.0021, Galaxies[q][0])
        plt.xlabel('xstd')
        plt.ylabel('ystd')
        plt.title('Standard Deviations of Centroid Calculations')
    pdf.savefig()
    plt.close()
