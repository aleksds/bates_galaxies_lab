#
# Author: Grayson Petter
# WISE.py
# Takes in a name of a galaxy and reads in the associated WISE fits file, calculates a flux in Jy and colors

import numpy as np
from astropy.table import Table
import os


#projpath = '/Users/graysonpetter/Desktop/IRSFRs/'
projpath = os.getcwd()+'/'

# convert from WISE magnitudes to fluxes in Jy
def mag_to_flux(name):

	WISEdir = projpath + 'unWISE/%s.fits' % name
	t = Table.read(WISEdir)

	if name == "J1506" or name == "J1219":
		print("Table", t)
	ra_col = t['ra']
	for i,r in enumerate(ra_col):
		if r == 226.65124:
			print ("Index of 1506:", i)

	# constants given at http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#example
	w_three_const = 31.674
	w_four_const = 8.363
	# calculate W3, W4 fluxes in Jy
	print("t[w3_mag]: ", t['w3_mag'])
	flux_three = w_three_const * (10 ** (-(float(t['w3_mag']) / 2.5)))
	flux_four = w_four_const * (10 ** (-(float(t['w4_mag']) / 2.5)))
	# errors computed by error propagation formula
	flux_three_err = w_three_const*np.log(10)/2.5*(10**(-(float(t['w3_mag']) / 2.5)))*float(t['w3_mag_err'])
	flux_four_err = w_four_const*np.log(10)/2.5*(10**(-(float(t['w4_mag']) / 2.5)))*float(t['w4_mag_err'])

	return flux_three, flux_four, flux_three_err, flux_four_err


# calculate WISE colors and errors
def colors(name):
	pat = projpath + 'unWISE/%s.fits' % name
	t = Table.read(pat)
	# W1-W2
	one_two = float(t['w1_mag'])-float(t['w2_mag'])
	one_two_err = np.sqrt((float(t['w1_mag_err']))**2+(float(t['w2_mag_err']))**2)
	# W3-W4
	three_four = float(t['w3_mag'])-float(t['w4_mag'])
	three_four_err = np.sqrt((float(t['w3_mag_err'])) ** 2 + (float(t['w4_mag_err'])) ** 2)
	# W2-W3
	two_three = float(t['w2_mag'])-float(t['w3_mag'])
	two_three_err = np.sqrt((float(t['w2_mag_err'])) ** 2 + (float(t['w3_mag_err'])) ** 2)
	return one_two, three_four, two_three, one_two_err, three_four_err, two_three_err

