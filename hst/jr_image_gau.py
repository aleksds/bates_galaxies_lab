# Josh Rines
# 20160923
#
# the goals of this code are to do the following:
# (1) read in three images taken in three different filters
# (2) plot a "postage stamp" version of each image, centered on the
#     science target
# (3) combine the three image to make a single RGB image (not yet implemented)
#
#     For future reference, here are some relevant links (from a
#     Google search for  "astropy make a three color image"):
## Josh Rines
# 20160923
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
import numpy as np
import os
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

# define the directory that contains the images
dir = '/Users/jrines/data/test/'

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

# create a PDF file for the plots    
with PdfPages('jr_image_gau.pdf') as pdf:

    fig = plt.figure()

    # read in the F475W image
    file = glob.glob(dir+'final_F4*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header

    # create a "postage stamp" image centered on the science target
    xcen = 3388.
    ycen = 3504.
    dx = 500
    dy = 500
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,1)
    plot_image()
    plt.title('F475W')

    # plot the gaussian

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
    #print('mean pixel value in the image: ', np.mean(pixels))
    #print('median pixel value in the image: ', np.median(pixels))
    
    # print the parameters of the best-fit gaussian function
    #print('background value from gaussian fit: {0:6.3f}'.format(params_fit[1]))
    #print('sigma from gaussian fit: {0:6.3f}'.format(params_fit[2]))
    
    # visualize the histogram and the fit
    #fig = plt.figure()
    ax = fig.add_subplot(2,3,4)
    ax.plot(bin_mids, hist, 'k:', color='blue', label='pixel histogram')
    ax.plot(bin_mids[use], hist[use], color='green', label='range for model fit')
    ax.plot(bin_mids, model, 'k--', color='red', label='Gaussian model')
    plt.tick_params(axis='both', which='major', labelsize=8)

    # added a legend, axis labels, title, text on plot
    #legend = ax.legend(fontsize=13)
    #plt.xlabel('Pixel values', fontsize=18)
    #plt.tick_params(axis='x', which='major', labelsize=16)
    #plt.xlim([-50., 50.])
    #plt.ylabel('Number of pixels', fontsize=18)
    #plt.tick_params(axis='y', which='major', labelsize=14)
    #plt.title(header['TARGNAME'])
    #str1 = str(r'$\mu=$')
    #str2 = str('{0:5.2f}'.format(params_fit[1]))
    #str3 = str(r', $\sigma=$')
    #str4 = str('{0:5.2f}'.format(params_fit[2]))
    #xloc = params_fit[1] - 4*params_fit[2]
    #yloc = params_fit[0]*0.8
    #ax.text(xloc, yloc, str1+str2+str3+str4, fontsize=14)
      
    # read in the F814W image
    file = glob.glob(dir+'final_F8*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    #xcen = 3629.
    #ycen = 4153.
    #dx = 500
    #dy = 500
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,2)
    plot_image()
    plt.title('F814W')

    # plot the gaussian

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
    #print('mean pixel value in the image: ', np.mean(pixels))
    #print('median pixel value in the image: ', np.median(pixels))
    
    # print the parameters of the best-fit gaussian function
    #print('background value from gaussian fit: {0:6.3f}'.format(params_fit[1]))
    #print('sigma from gaussian fit: {0:6.3f}'.format(params_fit[2]))
    
    # visualize the histogram and the fit
    #fig = plt.figure()
    ax = fig.add_subplot(2,3,5)
    ax.plot(bin_mids, hist, 'k:', color='blue', label='pixel histogram')
    ax.plot(bin_mids[use], hist[use], color='green', label='range for model fit')
    ax.plot(bin_mids, model, 'k--', color='red', label='Gaussian model')
    plt.tick_params(axis='both', which='major', labelsize=8)

    # added a legend, axis labels, title, text on plot
    #legend = ax.legend(fontsize=13)
    #plt.xlabel('Pixel values', fontsize=18)
    #plt.tick_params(axis='x', which='major', labelsize=16)
    #plt.xlim([-50., 50.])
    #plt.ylabel('Number of pixels', fontsize=18)
    #plt.tick_params(axis='y', which='major', labelsize=14)
    #plt.title(header['TARGNAME'])
    #str1 = str(r'$\mu=$')
    #str2 = str('{0:5.2f}'.format(params_fit[1]))
    #str3 = str(r', $\sigma=$')
    #str4 = str('{0:5.2f}'.format(params_fit[2]))
    #xloc = params_fit[1] - 4*params_fit[2]
    #yloc = params_fit[0]*0.8
    #ax.text(xloc, yloc, str1+str2+str3+str4, fontsize=14)

    # read in the F160W image
    file = glob.glob(dir+'final_F1*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
    # create a "postage stamp" image centered on the science target
    #xcen = 3629.
    #ycen = 4153.
    #dx = 500
    #dy = 500
    stamp = data[round(ycen-dy):round(ycen+dy), round(xcen-dx):round(xcen+dx)]

    # plot the "postage stamp"
    ax = fig.add_subplot(2,3,3)
    plot_image()
    plt.title('F160W')

    # plot the gaussian

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
    #print('mean pixel value in the image: ', np.mean(pixels))
    #print('median pixel value in the image: ', np.median(pixels))
    
    # print the parameters of the best-fit gaussian function
    #print('background value from gaussian fit: {0:6.3f}'.format(params_fit[1]))
    #print('sigma from gaussian fit: {0:6.3f}'.format(params_fit[2]))
    
    # visualize the histogram and the fit
    #fig = plt.figure()
    ax = fig.add_subplot(2,3,6)
    ax.plot(bin_mids, hist, 'k:', color='blue', label='pixel histogram')
    ax.plot(bin_mids[use], hist[use], color='green', label='range for model fit')
    ax.plot(bin_mids, model, 'k--', color='red', label='Gaussian model')
    plt.tick_params(axis='both', which='major', labelsize=8)

    # added a legend, axis labels, title, text on plot
    #legend = ax.legend(fontsize=13)
    #plt.xlabel('Pixel values', fontsize=18)
    #plt.tick_params(axis='x', which='major', labelsize=16)
    #plt.xlim([-50., 50.])
    #plt.ylabel('Number of pixels', fontsize=18)
    #plt.tick_params(axis='y', which='major', labelsize=14)
    #plt.title(header['TARGNAME'])
    #str1 = str(r'$\mu=$')
    #str2 = str('{0:5.2f}'.format(params_fit[1]))
    #str3 = str(r', $\sigma=$')
    #str4 = str('{0:5.2f}'.format(params_fit[2]))
    #xloc = params_fit[1] - 4*params_fit[2]
    #yloc = params_fit[0]*0.8
    #ax.text(xloc, yloc, str1+str2+str3+str4, fontsize=14)
    
    pdf.savefig()
    plt.close()
    os.system('open %s & ' % ' jr_image_gau.pdf')
