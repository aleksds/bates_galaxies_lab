#magnitudes to mass of stellar populations code
#10^(magnitude/(-2.5))*10^(9) = flux in nanomaggies
import numpy as np
from astropy.table import Table
import sys

#the code below allows one to provide an argument when using ipython so one can put in the desired folder (which is named based on the date and time)
date_time = sys.argv[1]
print(date_time)

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']



#SERSIC-INDEPENDENT
mags = np.zeros([12,2])
for w in range(0,12):
    for i in range(0,2):
        file = 'autogalfit/'+date_time+'_sersic_independent/'+galaxies[w]+'_'+filters[i]+'_coarse.galfit.01.band'
        with open(file) as f:
            content = f.readlines()
        mags[w][i] = np.float(content[47][4:10])

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(mags[w][0]/(-2.5))*10**9)
    iflux[w] = (10**(mags[w][1]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

table = Table([vflux, vivar, iflux, iivar], names=('sersic-indy 475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)



#SERSIC-SIMULTANEOUS
mags = np.zeros([12,2])
for w in range(0,12):
    file = 'autogalfit/'+date_time+'_sersic_simultaneous/'+galaxies[w]+'_F814W_F475W_sersic_output.galfit.01.band'
    with open(file) as f:
        content = f.readlines()
    mags[w][0] = np.float(content[47][11:17])
    mags[w][1] = np.float(content[47][4:10])

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(mags[w][0]/(-2.5))*10**9)
    iflux[w] = (10**(mags[w][1]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

table = Table([vflux, vivar, iflux, iivar], names=('sersic-simult 475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)



#PSF-INDEPENDENT
mags = np.zeros([12,2])
for w in range(0,12):
    for i in range(0,2):
        file = 'autogalfit/'+date_time+'_psf_independent/'+galaxies[w]+'_'+filters[i]+'_coarse.galfit.01.band'
        with open(file) as f:
            content = f.readlines()
        mags[w][i] = np.float(content[47][4:10])

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(mags[w][0]/(-2.5))*10**9)
    iflux[w] = (10**(mags[w][1]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

table = Table([vflux, vivar, iflux, iivar], names=('psf-indy 475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)



#PSF-SIMULTANEOUS
mags = np.zeros([12,2])
for w in range(0,12):
    file = 'autogalfit/'+date_time+'_psf_simultaneous/'+galaxies[w]+'_F814W_F475W_psf_output.galfit.01.band'
    with open(file) as f:
        content = f.readlines()
    mags[w][0] = np.float(content[47][11:17])
    mags[w][1] = np.float(content[47][4:10])

#flux in nanomaggies
vflux=np.zeros(12)
iflux=np.zeros(12)
for w in range(0,12):
    vflux[w] = (10**(mags[w][0]/(-2.5))*10**9)
    iflux[w] = (10**(mags[w][1]/(-2.5))*10**9)

#inverse variances
vivar=np.zeros(12)
iivar=np.zeros(12)
for i in range(0,12):
    vivar[i] = 1/(vflux[i]*0.05)**2
    iivar[i] = 1/(iflux[i]*0.05)**2

table = Table([vflux, vivar, iflux, iivar], names=('psf-simult 475 flux', '475 ivar', '814 flux', '814 ivar'))
print(table)
