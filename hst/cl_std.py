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
    return np.array([x1,x2,x3])
def plot_imageY(xcen,ycen):
    std = np.std(stamp[stamp==stamp])
    x1, y1 = centroid_com(stamp)
    x2, y2 = centroid_1dg(stamp)
    x3, y3 = centroid_2dg(stamp)
    return np.array([y1,y2,y3])
# define the directory that contains the images
dir = os.environ['HSTDIR']

# parameters for the initial guess at the galaxy centroid and the size
Galaxies = [['J0826', 3628.7, 4153.8], ['J0901', 3934, 4137], ['J0905', 3386.5, 3503.2], ['J0944', 3477.5, 3404.3], ['J1107', 3572.9, 3339.1], ['J1219', 3803., 4170.], ['J1341', 3885.6, 4164.3], ['J1506', 4149.2, 3921.7], ['J1558', 3787., 4186.], ['J1613', 4175., 3827.], ['J2116', 3566.9, 3435.9], ['J2140', 4067, 4054.4]]

dx = 5
dy = 5

xone = []
yone = []
xtwo = []
ytwo = []
xthree = []
ythree = []
xfour = []
yfour = []
xfive = []
yfive = []
xsix = []
ysix = []
xseven = []
yseven = []
xeight = []
yeight = []
xnine = []
ynine = []
xten = []
yten = []
xeleven = []
yeleven = []
xtwelve = []
ytwelve = []

# create a PDF file for the plots
yolo = ['centroid_com','centroid_1dg','centroid_2dg']
foureightone = [4,8,1]
collection = ['F475W','F814W','F160W']
colors = ['blue','green','red']

for j in range(0,len(Galaxies)):
        for i in range(0,len(collection)):
            file = glob.glob(dir+Galaxies[j][0]+'_final_F'+str(foureightone[i])+'*sci.fits')
            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
            stamp = data[round(Galaxies[j][2]-dy):round(Galaxies[j][2]+dy), round(Galaxies[j][1]-dx):round(Galaxies[j][1]+dx)]
        if j==0:
            xone = np.append(xone,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yone = np.append(yone,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==1:
            xtwo = np.append(xtwo,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            ytwo = np.append(ytwo,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==2:
            xthree = np.append(xthree,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            ythree = np.append(ythree,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==3:
            xfour = np.append(xfour,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yfour = np.append(yfour,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==4:
            xfive = np.append(xfive,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yfive = np.append(yfive,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==5:
            xsix = np.append(xsix,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            ysix = np.append(ysix,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==6:
            xseven = np.append(xseven,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yseven = np.append(yseven,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==7:
            xeight = np.append(xeight,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yeight = np.append(yeight,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==8:
            xnine = np.append(xnine,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            ynine = np.append(ynine,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==9:
            xten = np.append(xten,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yten = np.append(yten,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==10:
            xeleven = np.append(xeleven,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            yeleven = np.append(yeleven,plot_imageY(Galaxies[j][1],Galaxies[j][2]))
        if j==11:
            xtwelve = np.append(xtwelve,plot_imageX(Galaxies[j][1],Galaxies[j][2]))
            ytwelve = np.append(ytwelve,plot_imageY(Galaxies[j][1],Galaxies[j][2]))

stdxone = np.std(xone)
stdyone = np.std(yone)
stdxtwo = np.std(xtwo)
stdytwo = np.std(ytwo)
stdxthree = np.std(xthree)
stdythree = np.std(ythree)
stdxfour = np.std(xfour)
stdyfour = np.std(yfour)
stdxfive = np.std(xfive)
stdyfive = np.std(yfive)
stdxsix = np.std(xsix)
stdysix = np.std(ysix)
stdxseven = np.std(xseven)
stdyseven = np.std(yseven)
stdxeight = np.std(xeight)
stdyeight = np.std(yeight)
stdxnine = np.std(xnine)
stdynine = np.std(ynine)
stdxten = np.std(xten)
stdyten = np.std(yten)
stdxeleven = np.std(xeleven)
stdyeleven = np.std(yeleven)
stdxtwelve = np.std(xtwelve)
stdytwelve = np.std(ytwelve)

with PdfPages('cl_std.pdf') as pdf:
    fig = plt.figure()
    plt.scatter(stdxone,stdyone,label='J0826')
    plt.scatter(stdxtwo,stdytwo,label='J0901')
    plt.scatter(stdxthree,stdythree,label='J0905')
    plt.scatter(stdxfour,stdyfour,label='J0944')
    plt.scatter(stdxfive,stdyfive,label='J1107')
    plt.scatter(stdxsix,stdysix,label='J1219')
    plt.scatter(stdxseven,stdyseven,label='J1341')
    plt.scatter(stdxeight,stdyeight,label='J1506')
    plt.scatter(stdxnine,stdynine,label='J1558')
    plt.scatter(stdxten,stdyten,label='J1613')
    plt.scatter(stdxeleven,stdyeleven,label='J2116')
    plt.scatter(stdxtwelve,stdytwelve,label='J2140')
    plt.xlabel('x std',fontsize=14)
    plt.ylabel('y std',fontsize=14)
    plt.title('Standard Deviation of Centroid Calculations',fontsize=16)
    legend = plt.legend(loc='upper center',prop={'size':7})
    pdf.savefig()
    plt.close()
