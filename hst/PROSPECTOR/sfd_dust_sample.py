from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import numpy as np
from astropy.io import fits
import astropy.units as u
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


j2118dir = os.path.join(os.getenv('HIZEA_PROJECT'), 'j2118-nebula')

phot = fits.open(j2118dir+'/hizea_photo_galex_wise_v1.0.fit')

print(phot[1].data)

ra = phot[1].data['ra']
dec = phot[1].data['dec']

ebv = np.zeros(len(ra))
for i in range(0,len(ra)):
    coords = SkyCoord(ra[i], dec[i], frame='icrs', unit=u.deg)
    sfd = SFDQuery()
    ebv[i] = sfd(coords)
    print(phot[1].data['SHORT_NAME'][i], ebv[i], phot[1].data['EBV_SFD'][i])


filename = 'ebv_check.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    
    plt.scatter(ebv, phot[1].data['EBV_SFD'])
    plt.xlabel('E(B-V) from SFDQuery()')
    plt.ylabel('EBV_SFD from hizea_photo_galex_wise_v1.0.fit')

    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.scatter(ebv, (ebv-phot[1].data['EBV_SFD'])/ebv)
    plt.xlabel('E(B-V) from SFDQuery()')
    plt.ylabel('EBV_SFD % difference (SFDQuery - EBV_SFRD)/SFDQuery')
    plt.ylim([-0.05, 0.05])

    pdf.savefig()
    plt.close()


    fig = plt.figure()
    
    plt.scatter(ebv, (ebv-phot[1].data['EBV_MW'])/ebv)
    plt.xlabel('E(B-V) from SFDQuery()')
    plt.ylabel('EBV_MW % difference (SFDQuery - EBV_SFRD)/SFDQuery')
    #plt.ylim([-0.05, 0.05])

    pdf.savefig()
    plt.close()
    
os.system('open %s &' % filename)
