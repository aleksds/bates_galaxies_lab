# Aleks Diamond-Stanic
# ft. Sophia CW Gottlieb I
# 20160915
#
# written with the following goals:
# (1) quantify the distribution of pixel values in an image
# (2) estimate the mean and standard deviation of the sky background
# (3) fit a Gaussian function to the background pixel values

# import relevant Python modules
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

# read in an image file
dir = '/Users/sgottlie/data/test/'
file = glob.glob(dir+'final_F*sci.fits')
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
plt.plot(bin_mids, hist)
plt.plot(bin_mids[use], hist[use], color='green')
plt.plot(bin_mids, model, color='red')

plt.show()
