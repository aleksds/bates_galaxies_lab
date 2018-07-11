import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv

# Categories of Galaxies Outflows
#   1) clear outflow
#     1*) exciting outflow along major axis, inflow along minor axis
#   2) clear opposite outflow
#   3) no clear blue/red extraplanar separation
#
#   0) clear blue/red separation but no clear orientation

# read information from from database files
drp = fits.open('drpall-v2_4_3.fits')                                   # drpall sdss file mpl7 file
drpdata = drp[1].data

nsa = fits.open('1-nsa_v1_0_1.fits')                                    # nsa catalog
nsa_data = nsa[1].data

# dprall arrays
redshift = drpdata.field('z')                                           # redshift
ba = drpdata.field('nsa_sersic_ba')                                     # b/a ratio
stellar_mass = drpdata.field('nsa_sersic_mass')                         # stellar mass
absmag = drpdata.field('nsa_elpetro_absmag')                            # absolute magnitude
nuvi= absmag[:,1]-absmag[:,5]                                           # nuv-i
iabs = absmag[:,5]                                                      # i band absolute magnitude
plateifu = drpdata.field('plateifu')                                    # plate ifu array


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

# concentration value array
c_manga = th90_manga / th50_manga

# bool for location of MaNGA galaxies out of all galaxies
blah = drpdata.field('srvymode') == 'MaNGA dither'

# update concentration for just MaNGA galaxies
c_manga_gal = c_manga[blah]
ba_gal = ba[blah]

# define edge on and late-type galaxies within MaNGA set of galaxies
late = c_manga_gal < 2.6                                                # late-type concentration defined as less than 2.6
ba_gal_late = ba_gal[late]
edge = abs(ba_gal_late) < 0.3                                           # edge-on galaxies defined as b/a ration less than 0.3 radians
ba_gal_late_edge = ba_gal_late[edge]

# isolating data array of only MaNGA late-type-edge-on galaxies
good_plates = np.sort(plateifu[blah][late][edge])                       # Galaxies of interest
size = len(good_plates)                                                 # variable equal to the total number of MaNGA late-edge-on galaxies

# update concentration and ba arrays to include only late-type galaxies
stellar_mass_gal = np.sort(stellar_mass[blah][late][edge])
nuv_i = np.sort(nuvi[blah])
i_abs = np.sort(iabs[blah])
r_shift = np.sort(redshift[blah])
nuvi_gal = nuvi[blah][late][edge]
iabs_gal = iabs[blah][late][edge]
redshift_gal = redshift[blah][late][edge]


###### Importing data from csv file
# Empty lists for csv file data
outflow_assign = []                                                     # for the data on the file about outflow
br_assign = []                                                          # for the data on the file about b/r sep

# open file containing classification of late-type-edge-on galaxies
with open('em_outflow.csv', 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter=','):
        outflow_assign.append(row[3])
        br_assign.append(row[4])

del outflow_assign[0]                                                    # deleting the column label outflow the csv file data set
del br_assign[0]                                                     # deleting the column label blue/red sep from the csv file data set

# Empty lists for categorization
c1 = np.zeros(size)                                                           # clear outflow
c2 = np.zeros(size)                                                             # clear opposite outflow
c3 = np.zeros(size)                                                            # no clear blue/red separation
c0 = np.zeros(size)                                                         # clear blue/red but no clear orientation

# counters to see how many galaxies lie in each outflow category
g = 0                                                                  # category 1 counter
n = 0                                                                  # category 2 counter
ok = 0                                                                 # category 3 counter
uh = 0                                                                 # category 4 counter

# outflows category assignment
for i in range(0,size):
    if outflow_assign[i] == 'yes':
        c1[i] = 10

        g = g + 1
    elif outflow_assign[i] == 'no':
        c2[i] = 10
        n = n + 1
    elif outflow_assign[i] == 'unclear':
        c0[i] = 10
# blue/red separation category assignment
for b in range(0,size):
    if br_assign[b] == 'unclear':
        c3[b] = 10
        c0[b] = 0
        ok = ok + 1
    elif br_assign == 'yes ' and c0[b] == 10:
        c0[b] = 10

# outflow categories bool arrays
cat_1 = abs(np.array(c1)) > 0
cat_2 = abs(np.array(c2)) > 0
cat_3 = np.array(c3) > 0
cat_0 = np.array(c0) > 0

filename = 'outflow_nuvi_absmag_hist.pdf'
with PdfPages(filename) as pdf:

    # plot 1: All mpl-7 galaxies
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs[blah], nuvi[blah], color='lightsteelblue', label='All MPL-7 galaxies', s=5, alpha=1.0)
    plt.xlabel('i-band')
    plt.ylabel('NUV-i')
    plt.title('I-band Absolute Magnitude vs NUV-i of All MPL7 Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 2: all late type galaxies
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs[blah][late], nuvi[blah][late], color='cyan', label='Late-type Galaxies', s=5, marker='^', alpha=1.0)
    plt.xlabel('i-band')
    plt.ylabel('NUV-i')
    plt.title('I-band Absolute Magnitude vs NUV-i of Late-type Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 3: late-type-edge-on
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs, nuvi, color='navy', label='Late-type-edge-on Galaxies', s=5, alpha=1.0)
    plt.xlabel('i-band absolute magnitude')
    plt.ylabel('nuv-i')
    plt.title('I-band Absolute Magnitude vs NUV-i of Late-type-edge-on Galaxies')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 4: late-type-edge-on
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs_gal[cat_1], nuvi_gal[cat_1], color='red', label='Late-type-edge-on Galaxies', s=5, alpha=1.0)
    plt.xlabel('i-band absolute magnitude')
    plt.ylabel('nuv-i')
    plt.title('I-band Absolute Magnitude vs NUV-i of Late-type-edge-on Galaxies with clear outflows')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 5: composite
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs[blah], nuvi[blah], color='lightsteelblue', label='All MPL-7 galaxies', s=5, alpha=0.5)
    plt.scatter(iabs[blah][late], nuvi[blah][late], color='cyan', label='Late-type Galaxies', s=5, alpha=0.8)
    plt.scatter(iabs_gal, nuvi_gal, color='navy', label='Late-type-edge-on Galaxies', s=5, alpha=1.0)
    plt.xlabel('i-band absolute magnitude')
    plt.ylabel('nuv-i')
    plt.title('I-band Absolute Magnitude vs NUV-i')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 6: composite
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(1, 8)
    plt.scatter(iabs, nuvi, color='lightsteelblue', label='All MPL-7 galaxies', s=5, alpha=0.5)
    plt.scatter(iabs[blah][late], nuvi[blah][late], color='cyan', label='Late-type Galaxies', s=5, alpha=0.8)
    plt.scatter(iabs_gal, nuvi_gal, color='navy', label='Late-type-edge-on Galaxies', s=5, alpha=1.0)
    plt.scatter(iabs_gal[cat_1], nuvi_gal[cat_1], color='red', label='Clear outflows', s=5, alpha=1.0)
    plt.xlabel('i-band absolute magnitude')
    plt.ylabel('nuv-i')
    plt.title('I-band Absolute Magnitude vs NUV-i')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

os.system("open %s &" % filename)
