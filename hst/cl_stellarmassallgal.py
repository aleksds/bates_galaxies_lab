#magnitudes to mass of stellar populations code
#10^(magnitude/(-2.5))*10^(9) = flux in nanomaggies
import numpy as np
from astropy.table import Table
#J0826 magnitudes
vmag = [20.48,20.41,20.28,20.77,20.54,20.61,19.99,20.07,20.56,20.67,21,20.48]
imag = [20.55,20.43,20.62,20.65,20.56,20.16,20.31,20.39,20.27,20.07,20.76,20.82]

#J0826 flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(vmag[w]/(-2.5))*10**9)*2.40
    iflux[w] = (10**(imag[w]/(-2.5))*10**9)*4.11

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
