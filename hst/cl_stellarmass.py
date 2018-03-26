#magnitudes to mass of stellar populations code
#10^(magnitude/(-2.5))*10^(9) = flux in nanomaggies

from astropy.table import Table
#J0826 magnitudes
vmag = 20.48
imag = 20.55

#J0826 flux in nanomaggies
vflux = 10**(vmag/(-2.5))*10**9
iflux = 10**(imag/(-2.5))*10**9

#inverse variances
vivar = 1/(vflux*0.05)**2
iivar = 1/(iflux*0.05)**2

fsfflux = [vflux]
fsfivar = [vivar]
efflux = [iflux]
efivar = [iivar]
table = Table([fsfflux, fsfivar, efflux, efivar], names=('475 flux', '475 ivar', '814 flux', '814 ivar'))
