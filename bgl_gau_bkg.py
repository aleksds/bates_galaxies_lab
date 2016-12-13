# Aleks Diamond-Stanic
# 20160915 -- 20160921
#
# Quick description: This code analyzes the distribution of pixel
# values for an entire image and quantifies the mean background flux
# and its dispersion by fitting a Gaussian function to a binned
# histogram of pixel values.
#
# Current status: The current version analyzes the F814W image for the
# galaxy J0826.
#
# Future developments: Could expand to analyze all three filters for
# any individual galaxy.  Could loop over all galaxies.
#
# Original notes:
#
# written with the following goals:
# (1) quantify the distribution of pixel values in an image
# (2) estimate the mean and standard deviation of the sky background
# (3) fit a Gaussian function to the background pixel values

# import relevant Python modules
import numpy as np
from astropy.io import fits
import glob
import os
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin

# create a PDF file for the plots    
with PdfPages('bgl_gau_bkg.pdf') as pdf:

    # read in an image file
    dir = os.environ['HSTDIR']
    file = glob.glob(dir+'J0826*final_F814W*sci.fits')
    hdu = fits.open(file[0])
    data, header = hdu[0].data, hdu[0].header
    
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
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(bin_mids, hist, 'k:', color='blue', label='pixel histogram')
    ax.plot(bin_mids[use], hist[use], color='green', label='range for model fit')
    ax.plot(bin_mids, model, 'k--', color='red', label='Gaussian model')

    # added a legend, axis labels, title, text on plot
    legend = ax.legend(fontsize=13)
    plt.xlabel('Pixel values', fontsize=18)
    plt.tick_params(axis='x', which='major', labelsize=16)
    plt.xlim([-50., 50.])
    plt.ylabel('Number of pixels', fontsize=18)
    plt.tick_params(axis='y', which='major', labelsize=14)
    plt.title(header['TARGNAME'])
    str1 = str(r'$\mu=$')
    str2 = str('{0:5.2f}'.format(params_fit[1]))
    str3 = str(r', $\sigma=$')
    str4 = str('{0:5.2f}'.format(params_fit[2]))
    xloc = params_fit[1] - 4*params_fit[2]
    yloc = params_fit[0]*0.8
    ax.text(xloc, yloc, str1+str2+str3+str4, fontsize=14)
    
    pdf.savefig()
    plt.close()

    os.system('open %s &' % 'bgl_gau_bkg.pdf')
