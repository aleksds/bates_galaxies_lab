# Aleks Diamond-Stanic
# ft. Sophia CW Gottlieb I
# 20160915
#
# Sophia edits 20170627 for all galaxies, no more plotting
#
# written with the following goals:
# (1) quantify the distribution of pixel values in an image
# (2) estimate the mean and standard deviation of the sky background
# (3) fit a Gaussian function to the background pixel values

# import relevant Python modules
import os
import numpy as np
from astropy.io import fits
import glob
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# define a function for the gaussian model
def gaussian(x, peak, center, sigma):
    return peak * np.exp(((x-center)/sigma)**2*(-0.5))

# define a least squares objective function
def objective_function(params):
    peak, center, sigma = params
    residuals = (gaussian(bin_mids[use], peak, center, sigma) - hist[use]) / errors[use]
    return np.sum((residuals)**2) / nbin
w=0
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
# read in an image file
dir = os.environ['HSTDIR']#'/Users/sgottlie/data/test/'
sigmas = np.zeros(len(galaxies))

#def rSky(fileName):
for w in range(0,len(galaxies)):
    print(galaxies[w])
    file = glob.glob(dir+galaxies[w]+'_final_F*sci.fits')
    hdu = fits.open(file[2])
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
    
    sigmas[w]= params_fit[2]

