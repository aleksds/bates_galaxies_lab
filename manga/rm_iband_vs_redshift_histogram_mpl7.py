import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
 
drp = fits.open('drpall-v2_4_3.fits')
drpdata = drp[1].data

redshift = drpdata.field('z')
ba = drpdata.field('nsa_sersic_ba')
stellar_mass = drpdata.field('nsa_sersic_mass')
absmag = drpdata.field('nsa_elpetro_absmag')
nuv_i= absmag[:,1]-absmag[:,5]
iabs = absmag[:,5]
plateifu = drpdata.field('plateifu')

# read in information from NSA catalog
nsa = fits.open('1-nsa_v1_0_1.fits')
nsa_data = nsa[1].data


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

# get just the MaNGA galaxies
blah = drpdata.field('srvymode') == 'MaNGA dither'
print(len(c_manga[blah]))

# update the concentration and ba arrays
c_manga_gal = c_manga[blah]
ba_gal = ba[blah]
stellar_mass_gal = stellar_mass[blah]
nuv_i_gal = nuv_i[blah]
iabs_gal = iabs[blah]
redshift_gal = redshift[blah] 

# find late-type manga galaxies
late = c_manga_gal < 2.6

# update concentration and ba arrays
ba_gal_late = ba_gal[late]
c_manga_gal_late = c_manga_gal[late]
stellar_mass_gal_late = stellar_mass_gal[late]
nuv_i_gal_late = nuv_i_gal[late]
iabs_gal_late = iabs_gal[late]
redshift_gal_late = redshift_gal[late]

# find edge-one late-type galaxies
edge = ba_gal_late < 0.3

# update concentration and ba arrays
c_manga_gal_late_edge = c_manga_gal_late[edge]
ba_gal_late_edge = ba_gal_late[edge]
stellar_mass_gal_late_edge = stellar_mass_gal_late[edge]
nuv_i_gal_late_edge = nuv_i_gal_late[edge]
iabs_gal_late_edge = iabs_gal_late[edge]
redshift_gal_late_edge = redshift_gal_late[edge]



filename = 'iband_v_redshift_histogram.pdf'
with PdfPages(filename) as pdf:

    #plot 1: b/a ratio

    #log_stellar_mass_gal = np.log10(stellar_mass_gal)
    #use = np.isfinite(log_stellar_mass_gal)

    #log_stellar_mass_gal_late_edge = np.log10(stellar_mass_gal_late_edge)
    #yo = np.isfinite(log_stellar_mass_gal_late_edge)

    #log_stellar_mass_gal_late = np.log10(stellar_mass_gal_late)
    #mo = np.isfinite(log_stellar_mass_gal_late)

    
    #plot 1: All MPL-7
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.16)
    plt.scatter(iabs_gal, redshift_gal, color='lightsteelblue', label= 'All MPL-7 galaxies', s=5, alpha = 1.0)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('I-band Absolute Magnitude vs Redshift of All MPL7 Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    
    #plot 2: all late type galaxies
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.16)
    plt.scatter(iabs_gal_late, redshift_gal_late, color='cyan', label= 'Late-type Galaxies', s=5, marker= '^', alpha = 1.0)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('I-band Absolute Magnitude vs Redshift of Late-type Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    #plot 3: late-type-edge-on
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.16)
    plt.scatter(iabs_gal_late_edge, redshift_gal_late_edge, color='navy', label= 'Late-type-edge-on Galaxies', s=5, alpha = 1.0)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('I-band Absolute Magnitude vs Redshift of Late-type-edge-on Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    #plot 4: composite
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.16)
    plt.scatter(iabs_gal, redshift_gal, color='lightsteelblue', label= 'All MPL-7 galaxies', s=5, alpha = 0.5)
    plt.scatter(iabs_gal_late, redshift_gal_late, color='cyan', label= 'Late-type Galaxies', s=5, alpha = 0.8)
    plt.scatter(iabs_gal_late_edge, redshift_gal_late_edge, color='navy', label= 'Late-type-edge-on Galaxies', s=5, alpha = 1.0)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close


    
    
os.system("open %s &" % filename)
    
