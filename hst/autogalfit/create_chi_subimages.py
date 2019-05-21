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
fits_files = glob.glob(dir+'/'+gal+'*'+band+'*output.fits')
band_files = glob.glob(dir+'/'+gal+'*'+band+'*band')
subimages_chivalues = open('Chivalues_of_subimages.dat', 'a+')
subimages_chivalues.write("Subimages_chi_values\n")
mag = np.zeros(len(band_files))
sizepix = np.zeros(len(band_files))
chi_from_our_calculations = np.zeros(len(band_files))

xcen = 100
ycen = xcen
dx = 10
dy = dx

for i in range(0, len(fits_files)):
    #print(fits_files[i], band_files[i])
    hdu = fits.open(fits_files[i])
    data, header = hdu[0].data, hdu[0].header
    with open(band_files[i]) as f:
        content = f.readlines()
        unc, unc_head = hdu[4].data, hdu[4].header
        res, res_head = hdu[2].data, hdu[2].header

        unc_image_stamp = unc[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]
        res_image_stamp = res[(ycen-dy):(ycen+dy), (xcen-dx):(xcen+dx)]

        #print(res_image_stamp)
        #print(unc_image_stamp)
        #chi_from_our_calculations_for_v[i] = np.sum((res_image_stamp/unc_image_stamp)**2)/((101**2)*3)
        #chi_from_our_calculations_for_u[i] = np.sum((res_image_stamp/unc_image_stamp)**2)/((101**2)*3)
        chi_from_our_calculations[i] = np.sum((res_image_stamp/unc_image_stamp)**2)/((len(res_image_stamp)**2)*1)
        mag[i] = np.float(content[47][4:13])
        sizepix[i] = np.float(content[48][4:13])
        #print(chi_from_our_calculations[i])
        #print(chi_from_our_calculations_for_u[i])

mag_1d = np.unique(mag)
sizepix_1d = np.unique(sizepix)
chi_2d = np.zeros([len(mag_1d), len(sizepix_1d)])
for i in range(0,len(mag_1d)):
    for j in range(0, len(sizepix_1d)):
        test = np.where((mag == mag_1d[i]) & (sizepix == sizepix_1d[j]))
        print(mag[test], sizepix[test], chi_from_our_calculations[test])
        chi_2d[i,j] = chi_from_our_calculations[test[0][0]]

name = 'chi_values_'+gal+'_'+band+'.pdf'
file = 'chi_values_'+gal+'_'+band+'.dat'
file_object = open('chi_values_'+gal+'_'+band+'.dat','a+')
file_object.write('Chi Value from Galfit'+gal+'_'+band+'\n')
with PdfPages(name) as pdf:
    fig = plt.figure()

    plt.scatter(mag, sizepix, marker='o', color='orange')
    plt.yscale('log')
    plt.ylim([0.1, 3.])
    plt.xlabel('Half-Light radius in pixels')
    plt.ylabel('Magnitude')
    plt.title(dir+'/'+gal+'_'+band)

    for i in range(0, len(band_files)):
        plt.text(mag[i], sizepix[i], str(chi_from_our_calculations[i]), fontsize=5)

    pdf.savefig()
    plt.close()

    fig = plt.figure()
    print(mag_1d, sizepix_1d, chi_2d)
    plt.contourf(sizepix_1d, mag_1d, chi_2d, 20, cmap='RdGy')
    plt.colorbar()

    plt.xlabel('Half-Light radius in pixels')
    plt.ylabel('Magnitude')
    pdf.savefig()
    plt.close()

    fig = plt.figure()

    chi_min = np.min(chi_2d)
    #levels = np.array([1.0, 4.00, 9.00])+chi_min
    levels = np.array([2.3, 4.61, 9.21])+chi_min
    cs = plt.contour(sizepix_1d, mag_1d, chi_2d, levels, colors=['blue', 'green', 'red'])
    plt.clabel(cs, inline=1, fontsize=14)
    plt.xlim([0,np.max(sizepix)*0.5])
    #plt.ylim([0, 2])
    plt.xlabel('Half-Light radius in pixels', fontsize=14)
    plt.ylabel('Magnitude', fontsize=14)
    labels = ['68%', '90%', '99%']
    #for i in range(len(labels)):
    #    cs.collections[i].set_label(labels[i])

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
    #Prints value of the min and max half life radii
    print(np.min(x), np.max(x))
    file_object.write("Min x:"+str(np.min(x))+" Max x:"+str(np.max(x))+"\n")
    #Prints the average value and then the uncertainty value
    print((np.min(x) + np.max(x))/2, (np.max(x) - np.min(x))/2)
    file_object.write("Average Value of x:"+str(np.max(x)/2)+'\n' "Uncertainty Value:"+str(np.max(x) - np.min(x)/2)+'\n')
    #Prints value of the min and max mag values
    print(np.min(y), np.max(y))
    file_object.write('Min y:'+str(np.min(y))+' Max y:'+str(np.max(y))+'\n')
    #Prints the average value and then the uncertainty value 
    print((np.min(y) + np.max(y))/2, (np.max(y) - np.min(y))/2)
    file_object.write('Average Value of y:'+str(np.min(y)+np.max(y)/2)+'\n' 'Uncertainty Value:'+str(np.max(y) - np.min(y)/2)+'\n')
    #Prints the converted pixels to arcesc with the knowledge of arcsecs of 
    #of the sky on Earth based on redshift
    print(np.min(x)*0.025, np.max(x)*0.025)
    file_object.write('Arcsecs for min x:'+str(np.min(x)*0.025)+'\n' 'Arcsecs for max x:'+str(np.max(x)*0.025))
    plt.legend(loc='upper right', prop={'size': 15})

    pdf.savefig()
    plt.close()

os.system('open %s &' % name)
