#trying to extract wise match information specifically their magnitude and uncertainties in the infrared

from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
#data = np.zeros()

#extracting data from fits files
for i in range(0, len(galaxies)):

    files = glob.glob('J*.fits')
    print(galaxies[i])
    hdu = fits.open(files[i])
    print(hdu)

    data = hdu[1].data
    print('data:',data)

    w1_mag = hdu[1].data.field('w1_mag')
    print('w1_mag',w1_mag)
    w1_mag_err = hdu[1].data.field('w1_mag_err')
    print('w1_mag_err',w1_mag_err)

    w2_mag = hdu[1].data.field('w2_mag')
    print('w2_mag', w2_mag)
    w2_mag_err = hdu[1].data.field('w2_mag_err')
    print('w2_mag_err',w2_mag_err)

    w3_mag = hdu[1].data.field('w3_mag')
    print('w3_mag', w3_mag)
    w3_mag_err = hdu[1].data.field('w3_mag_err')
    print('w3_mag_err',w2_mag_err)

    w4_mag = hdu[1].data.field('w4_mag')
    print('w4_mag',w4_mag)
    w4_mag_err = hdu[1].data.field('w4_mag_err')
    print('w4_mag_err',w4_mag_err)
