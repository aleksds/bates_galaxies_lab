#magnitudes to mass of stellar populations code
#10^(magnitude/(-2.5))*10^(9) = flux in nanomaggies
import numpy as np
from astropy.table import Table
#magnitudes
vmag = [19.494,19.468,19.54,19.819,19.394,19.389,19.048,19.112,19.498,19.719,19.905,20.329]
imag = [18.647,18.496,18.984,19.038,18.569,18.155,18.544,18.642,18.143,18.5,19.104,18.934]

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(vmag[w]/(-2.5))*10**9)
    iflux[w] = (10**(imag[w]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

#fsfflux = [vflux]
#fsfivar = [vivar]
#efflux = [iflux]
#efivar = [iivar]
table = Table([vflux, vivar, iflux, iivar], names=('475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)




#now for psf fits
vmag = [19.686,19.627,19.647,19.929,19.699,19.995,19.192,19.183,19.949,20.598,20.262,20.431]
imag = [19.253,19.138,19.19,19.263,19.427,19.137,18.96,19.048,19.661,19.63,19.601,19.63]

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(vmag[w]/(-2.5))*10**9)
    iflux[w] = (10**(imag[w]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

#fsfflux = [vflux]
#fsfivar = [vivar]
#efflux = [iflux]
#efivar = [iivar]
table = Table([vflux, vivar, iflux, iivar], names=('475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)
