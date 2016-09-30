# Sophia CW Qottlieb I
# 20160923
#
# MODIFIED FROM IMSTMP
#
# the goals of this code are to do the following:
# (1) read in three images taken in three different filters
# (2) plot a "postage stamp" version of each image, centered on the
#     science target
# (3) combine the three image to make a single RGB image (not yet implemented)
#
#     For future reference, here are some relevant links (from a
#     Google search for  "astropy make a three color image"):
#
#     (A) "Making RGB images from FITS files with python / matplotlib"
#     http://www.astrobetter.com/blog/2010/10/22/making-rgb-images-from-fits-files-with-pythonmatplotlib/
#
#     (B) "Making a publication quality image"
#     https://python4astronomers.github.io/intro/quick-tour.html

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import minimize

# define a function to plot "postage stamp" images
def plot_image():
    std = np.std(stamp[stamp==stamp])
    plt.imshow(stamp, interpolation='nearest', origin = 'lower', vmin = -1.*std, vmax = 3.*std, cmap='bone')
    plt.tick_params(axis='both', which='major', labelsize=8)
    
# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

# define the directory that contains the images
dir = '/Users/sgottlie/data/test/'
filename = 'sg_gau_stamp.pdf'

# create a PDF file for the plots    
with PdfPages(filename) as pdf:
    

    fig = plt.figure()

    y = -1
    collection = ['F475W','F814W','F160W']
    for x in collection:
    # for x in range(0, len(collection)):

        y += 1
        
        # read in image
        # how can I modify this code to work for a loop?
        file = glob.glob(dir+'final_' + collection[y] + '*sci.fits')
        print(file)
        hdu = fits.open(file[0])
        data, header = hdu[0].data, hdu[0].header

        # create postage stamp image centered on target
        xcen = 3629.
        ycen = 4153.
        dx = 500
        dy = 500
        stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

        # plot the "postage stamp"
        # fig = plt.figure()
        ax = fig.add_subplot(2,3,y+1)
        plot_image()
        plt.title(x)

        # select all pixels with meaningful information and put them into a
        # one-dimensional array
        pixels = data[(data==data)].ravel()

        # generate a histogram of pixel values for this array
        bin_lo = -100
        bin_hi = 100
        bin_edges = np.linspace(bin_lo, bin_hi, 5000)
        hist, edges = np.histogram(pixels, bins=bin_edges)
        bin_mids = (bin_edges[0:-1] +  bin_edges[1:]) / 2.
        
        # decide on range of pixel values for gaussian fit
        rlo = -5
        rhi = 5
        use = (bin_mids > rlo) & (bin_mids < rhi)
        nbin = len(use)
        
        # estimate the errors in each bin by assuming poisson statistics: sqrt(n)
        errors = np.sqrt(hist)
        
        # decide on initial parameters for the gaussian
        params_initial = [1e4, 0, 10]
        
        # find the best-fit gaussian
        results = minimize(objective_function, params_initial, method='Powell')
        params_fit = results.x
        model = gaussian(bin_mids, *params_fit)
        
        # print the mean and median values for pixels in the image
        print('mean pixel value in the image: ', np.mean(pixels))
        print('median pixel value in the image: ', np.median(pixels))
        
        # print the parameters of the best-fit gaussian function
        print('background value from gaussian fit: {0:6.3f}'.format(params_fit[1]))
        print('sigma from gaussian fit: {0:6.3f}'.format(params_fit[2]))
        
        # visualize the histogram and the fit
        plt.plot(bin_mids, hist)
        plt.plot(bin_mids[use], hist[use], color='green')
        plt.plot(bin_mids, model, color='red')
        ax = fig.add_subplot(2,3,y+4)
        plot_image()
        plt.title(x)
        # plt.show()
        

        pdf.savefig()
        plt.close()
os.system('open %s &' % filename)

    ####### stuff and things
    # Aleks Diamond-Stanic
# ft. Sophia CW Gottlieb I
# 20160915
#
# written with the following goals:
# (1) quantify the distribution of pixel values in an image
# (2) estimate the mean and standard deviation of the sky background
# (3) fit a Gaussian function to the background pixel values

# import relevant Python modules


#
## read in an image file
## dir = '/Users/adiamond/data/20150203_HST/J0826+4305/coarse/F814W/'
## file = glob.glob(dir+'final*sci.fits')
## dir = '/Users/sgottlie/data/test/'
#file = glob.glob(dir+'final_F*sci.fits')
#hdu = fits.open(file[0])
#data, header = hdu[0].data, hdu[0].header
#
## select all pixels with meaningful information and put them into a
## one-dimensional array
#pixels = data[(data==data)].ravel()
#
## generate a histogram of pixel values for this array
#bin_lo = -100
#bin_hi = 100
#bin_edges = np.linspace(bin_lo, bin_hi, 5000)
#hist, edges = np.histogram(pixels, bins=bin_edges)
#bin_mids = (bin_edges[0:-1] +  bin_edges[1:]) / 2.
#
## decide on range of pixel values for gaussian fit
#rlo = -5
#rhi = 5
#use = (bin_mids > rlo) & (bin_mids < rhi)
#nbin = len(use)
#
## estimate the errors in each bin by assuming poisson statistics: sqrt(n)
#errors = np.sqrt(hist)
#
## decide on initial parameters for the gaussian
#params_initial = [1e4, 0, 10]
#
## find the best-fit gaussian
#results = minimize(objective_function, params_initial, method='Powell')
#params_fit = results.x
#model = gaussian(bin_mids, *params_fit)
#
## print the mean and median values for pixels in the image
#print('mean pixel value in the image: ', np.mean(pixels))
#print('median pixel value in the image: ', np.median(pixels))
#
## print the parameters of the best-fit gaussian function
#print('background value from gaussian fit: {0:6.3f}'.format(params_fit[1]))
#print('sigma from gaussian fit: {0:6.3f}'.format(params_fit[2]))
#
## visualize the histogram and the fit
#plt.plot(bin_mids, hist)
#plt.plot(bin_mids[use], hist[use], color='green')
#plt.plot(bin_mids, model, color='red')
#ax = fig.add_subplot(2,3,y+1)
#plot_image()
#plt.title(x)
#plt.show()
#    
