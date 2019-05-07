# Purpose: Is to determine the chi values from a subimage of the
# images given. We do this by taking the subimage of a galaxy
# starting at the center and moving out
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
gal = sys.argv[2]
band = sys.argv[3]
files = glob.glob(dir+'/*output.txt')
subimages_chivalues = open('Chivalues_of_subimages.dat', 'a+')
subimages_chivalues.write("Subimages_chi_values\n")
mag = np.zeros(len(files))
sizepix = np.zeros(len(files))

xcen = 50
ycen = 50
dx = 5
dy = 5

for i in range(0, len(files)):
    hdu = fits.open(files[0])
    data, header = hdu[0].data, hdu[0].header
    with open(files[i]) as f:
        # F814W corresponds roughly to rest-frame U-band
        vunc, vunc_head = hdu[12].data, hdu[12].header
        vres, vres_head = hdu[6].data, hdu[6].header
        # F475W corresponds roughly to rest-frame U-band
        uunc, uunc_head = hdu[13].data, hdu[13].header
        ures, ures_head = hdu[8].data, hdu[8].header

        vunc_image_stamp = vunc[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]
        vres_image_stamp = vres[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]
        uunc_image_stamp = uunc[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]
        ures_image_stamp = ures[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]
        chi_from_our_calculations_for_v[i] = np.sum((vres_image_stamp/vunc_image_stamp)**2)/((101**2)*3)
        chi_from_our_calculations_for_u[i] = np.sum((ures_image_stamp/uunc_image_stamp)**2)/((101**2)*3)
        mag[i] = np.float(content[47][4:13])
        sizepix[i] = np.float(content[48][4:13])
        print(chi_from_our_calculations_for_v[i])
        print(chi_from_our_calculations_for_u[i])

mag_1d = np.unique(mag)
sizepix_1d = np.unique(sizepix)
chi_2d = np.zeroes([len(mag_1d), len(size_1d)])
for i in range(0,len(mag_1d)):
    for j in range(0, len(sizepix_1d)):
        test = np.where((mag == mag_1d[i]) & (sizepix == sizepix_1d[j]))
        print(mag[test], sizepix[test], chi[test])
        chi_2d[j,i] = chi_from_our_calculations_for_v[test[0]]

name = 'chi_values_'+gal+'_'+band+'.pdf'

with PdfPages(name) as pdf:
    fig = plt.figure()

    plt.scatter(mag, sizepix, marker='o', color='orange')
    plt.yscale('log')
    plt.ylim([0.1, 3.])
    plt.ylabel('Half-Light radius in pixels')
    plt.xlabel('Magnitude')
    plt.title(dir+'/'+gal+'_'+band)

    for i in range(0, len(files)):
        plt.text(mag[i], sizepix[i], str(chi_from_our_calculations_for_v[i]), fontsize=5)

    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.contourf(mag_1d, sizepix_1d, chi_2d, 20, cmap='RdGy')
    plt.colorbar()

    plt.ylabel('Half-Light radius in pixels')
    plt.xlabel('Magnitude')
    pdf.savefig()
    plt.close()

    fig = plt.figure()

    chi_min = np.min(chi_2d)
    levels = np.array([1.0, 4.00, 9.00])+chi_min
    cs = plt.contour(mag_1d, sizepix_1d, chi_2d, levels, colors=['blue', 'green', 'red'])
    plt.clabel(cs, inline=1, fontsize=14)
    plt.ylim([0, 2])
    plt.ylabel('Half-Light radius in pixels', fontsize=14)
    plt.xlabel('Magnitude', fontsize=14)
    labels = ['68%', '95%', '99.7%']
    for i in range(len(labels)):
        cs.collections[i].set_label(labels[i])

    pzero = cs.collections[0].get_paths()[0]
    vzero = pzero.vertices
    xzero = vzero[:,0]
    yzero = vzero[:,1]

    if (len(cs.collections[0].get_paths()) > 1):
        pone = cs.collections[0].get_paths()[1]
        vone = pone.vertices
        xone = vone[:,0]
        yone = vone[:,1]

        x = np.append(xzero, xone)
        y = np.append(yzero, yone)
    else:
        x = xzero
        y = yzero
    print(x,y)
    print(np.min(x), np.max(x))
    print((np.min(x) + np.max(x))/2, (np.max(x) - np.min(x))/2)
    print(np.min(y), np.max(y))
    print((np.min(y) + np.max(y))/2, (np.max(y) - np.min(y))/2)
    print(np.min(y)*0.025, np.max(y)*0.025)
    plt.legend(loc='upper right', prop={'size': 15})

    pdf.savefig()
    plt.close()

os.system('open %s &' % name)
