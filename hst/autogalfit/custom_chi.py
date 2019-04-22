# Cristopher Thompson
# Date: 2019/04/22
# This code was creating to determine the accurate chi values for each of the galaxies.

import numpy as np
import sys
import os
import glob
from astropy.io import fits

dir = sys.argv[1]
files = glob.glob(dir+'/*output.fits')

for i in range (0, len(files)):
   raw_data=fits.open(files[i])
   vdat, vadt_head = raw_data[0].data, raw_data[0],header
   vunc, vunc_head = raw_data[12].data, raw_data[12].header
   vmod, vmod_head = raw_data[3].data, raw_data[3].header
   vres, vres_head = raw_data[6].data, raw_data[6].header
   n = np.sum((vres/vunc)**2)/(101**2))
   chi = n
   print(chi)

os.system('open %s' %name)
