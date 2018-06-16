import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
 
drp = fits.open('drpall-v2_4_3.fits')
drpdata = drp[1].data 


ba = drpdata.field('nsa_sersic_ba')
stellar_mass = drpdata.field('nsa_sersic_mass')
absmag = drpdata.field('nsa_elpetro_absmag')
gr= absmag[:,3] - absmag[:,4]
filename = 'histograms.pdf'

# read in information from NSA catalog
nsa = fits.open('1-nsa_v1_0_1.fits')
nsa_data = nsa[1].data

# check on a galaxy of interest
plate = 7977
galaxy = 12704
match = np.where((drpdata.field('plate') == plate) & (drpdata.field('ifudsgn') == str(galaxy)))
mangaid = drpdata[match].field('mangaid')[0]
nsaid = drpdata[match].field('nsa_nsaid')[0]

yo = np.where(nsa_data.field('NSAID') == nsaid)
name = nsa_data[yo[0][0]].field('IAUNAME')
print(name)
blah = drpdata[match].field('nsa_iauname')[0]
print(blah)

# match the two catalogs on the basis of right ascension and declination
manga_ra = drpdata.field('objra')
manga_dec = drpdata.field('objdec')
nsa_ra = nsa_data.field('RA')
nsa_dec = nsa_data.field('DEC')
c = SkyCoord(ra=manga_ra*u.degree, dec=manga_dec*u.degree)
catalog = SkyCoord(ra=nsa_ra*u.degree, dec=nsa_dec*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

# define th50, th90, and concentration arrays
th50 = nsa_data.field('ELPETRO_TH50_R')
th90 = nsa_data.field('ELPETRO_TH90_R')

th50_manga = th50[idx]
th90_manga = th90[idx]
c_manga = th90_manga / th50_manga

with PdfPages(filename) as pdf:

    #plot 1: b/a ratio
    fig = plt.figure()
    good = ba > 0
    bad_ba = ba > 0.3
    print(np.count_nonzero(bad_ba))
    print(np.shape(drpdata))
    plt.hist(bad_ba[good], color='teal', bins=np.arange(0.1,1.02,.02))
    plt.xlabel('b/a ratio')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of b/a Values in MaNGA sample')
    plt.axvline(x=0.3, color='k', linestyle='dashed')

    pdf.savefig()
    plt.close

    #plot 2: stellar mass
    fig = plt.figure()
    print(np.shape(drpdata))
    plt.hist(np.log10(stellar_mass[good]),color='teal', bins=np.arange(8.0,12.0,.2))
    plt.xlim(8,12)
    plt.xlabel('Stellar Mass (logM☉)')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of Stellar Mass in MaNGA sample')

    pdf.savefig()
    plt.close

    #plot 3: g-r
    fig = plt.figure()
    print(np.shape(drpdata))
    plt.hist(gr[good], bins=np.arange(0.0,1.0,0.05), color='teal')
    plt.xlabel('g-r')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of g-r in MaNGA sample')

    pdf.savefig()
    plt.close

    #concentration
    bad_c = c_manga > 2.6
    fig = plt.figure()
    print(np.shape(drpdata))
    plt.hist(bad_c[good], bins=np.arange(1.5,4.1,0.1), color='teal')
    plt.xlabel('Concentration')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of Concentration in MaNGA sample')

    pdf.savefig()
    plt.close
    
    
os.system("open %s &" % filename)
    
