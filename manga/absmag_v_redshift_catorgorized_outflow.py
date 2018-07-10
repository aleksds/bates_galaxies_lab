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

mpl7_dir = os.environ['MANGADIR_MPL7']
drp = fits.open(mpl7_dir + 'drpall-v2_4_3.fits')                        # drpall sdss file mpl7 file
drpdata = drp[1].data

nsa = fits.open(mpl7_dir + '1-nsa_v1_0_1.fits')                         # nsa catalog
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
nuv_i_gal = nuvi[blah][late][edge]
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
c1 = np.zeros(size)                                                             # clear outflow
c2 = np.zeros(size)                                                             # clear opposite outflow
c3 = np.zeros(size)                                                             # no clear blue/red separation
c0 = np.zeros(size)                                                             # clear blue/red but no clear orientation

# counters to see how many galaxies lie in each outflow category
g = 0                                                                  # category 1 counter
n = 0                                                                  # category 2 counter
ok = 0                                                                 # category 3 counter
uh = 0                                                                 # category 4 counter

# outflows category assignment
for b in range(0,size):
    if outflow_assign[b] == 'yes':
        c1[b] = 10
        g = g + 1
    elif outflow_assign[b] == 'no':
        c2[b] = 10
        n = n + 1
    elif outflow_assign[b] == 'unclear':
        c0[b] = 10
    
# blue/red separation category assignment
for b in range(0,size):
    if br_assign[b] == 'unclear':
        c3[b] = 10
        c0[b] = 10
        ok = ok + 1
    elif br_assign == 'yes ' and c0[b] == 10:
        c4[b] = 10

# outflow categories bool arrays
cat_1 = abs(np.array(c1)) > 0
cat_2 = abs(np.array(c2)) > 0
cat_3 = np.array(c3) > 0
cat_0 = np.array(c0) > 0


# PDF and Graphs

filename = 'outflow_hist.pdf'
with PdfPages(filename) as pdf:
# Abs mag v. Redshift
    # plot 1: all MPL-7 Late Type Edge on Galaxies
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.14)
    plt.scatter(iabs_gal, redshift_gal, color='teal', label= 'All ', s=5, alpha = 1.0)
    plt.xlabel('i-band absolute magnitude')
    plt.ylabel('redshift')
    plt.title('All MPL-7 Late-type-Edge-on Galaxies')
    plt.suptitle('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper right', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 2: good outflows -- category 1
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.14)
    plt.scatter(iabs_gal[cat_1], redshift_gal[cat_1], color='crimson', marker='*', label= 'galaxies w/ clear outflow', s=5, alpha = .8)
    plt.xlabel('i-Band Absolute Magnitude')
    plt.ylabel('Redshift')
    plt.title('Clear Outflow')
    plt.suptitle('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper right', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # plot 3: unclear outflows -- category 2
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0, 0.14)
    plt.scatter(iabs_gal[cat_2], redshift_gal[cat_2], color='black', marker='^', label='galaxies w/ clear opposite outflows', s=5, alpha = .8)
    plt.xlabel('i-Mand Absolute Magnitude')
    plt.ylabel('Redshift')
    plt.title('Clear Opposite Outflow Patterns')
    plt.suptitle('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper right', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close


    # plot 4: no outflows -- categories 0 and 3
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0, 0.14)
    plt.scatter(iabs_gal[cat_3], redshift_gal[cat_3], color='darkblue', label= 'Galaxies with no clear blue/red sep', s=5, alpha = .8)
    plt.scatter(iabs_gal[cat_0], redshift_gal[cat_0], color='springgreen', label='Galaxies with clear b/r but unclear outflow', s=5, alpha=.8)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('Unclear Outflow Galaxies')
    plt.suptitle('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper right', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close



    # plot 5: composite of all categories
    fig = plt.figure()
    plt.xlim(-16, -24)
    plt.ylim(0,0.14)
    plt.scatter(iabs_gal[cat_3], redshift_gal[cat_3], color='darkblue', label= 'Galaxies with no clear blue/red sep', s=5, alpha = .5)
    plt.scatter(iabs_gal[cat_0], redshift_gal[cat_0], color='springgreen', label='Galaxies with clear b/r but unclear outflow', s=5, alpha=.5)
    plt.scatter(iabs_gal[cat_2], redshift_gal[cat_2], color='black', marker='^', label='galaxies w/ clear opposite outflows', s=5, alpha=.8)
    plt.scatter(iabs_gal[cat_1], redshift_gal[cat_1], color='crimson', marker='*', label= 'galaxies w/ clear outflow', s=5, alpha = .8)
    plt.xlabel('i-band Absolute Magnatitude')
    plt.ylabel('Redshift')
    plt.title('I-band Absolute Magnitude vs Redshift')
    plt.legend(loc='upper right', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close


#     #plot 5: All MPL-7
#     fig = plt.figure()
#     plt.xlim(-16, -24)
#     plt.ylim(1,8)
#     plt.scatter(iabs_gal, nuv_i_gal, color='lightsteelblue', label= 'All MPL-7 MaNGA Ltyepe Eon galaxies', s=5, alpha = 1.0)
#     plt.xlabel('i-band')
#     plt.ylabel('NUV-i')
#     plt.title('I-band Absolute Magnitude vs NUV-i of All MaNGA Ltype Eon Galaxies')
#     plt.legend(loc='upper left', frameon=False, fontsize='x-small')
#
#     pdf.savefig()
#     plt.close
#
#     # Plot 6
#     fig = plt.figure()
#     plt.xlim(-16, -24)
#     plt.ylim(1, 8)
#     plt.scatter(iabs_gal[clear], nuv_i_gal[clear], color='red', marker='^',  label='Clear outflow', s=3 ,alpha=1.0)
#     plt.xlabel('i-band')
#     plt.ylabel('NUV-i')
#     plt.title('I-band Absolute Magnitude vs NUV-i of Galaxies with Clear Outflow Patterns')
#     plt.legend(loc='upper left', frameon=False, fontsize='x-small')
#
#     pdf.savefig()
#     plt.close
#
#     # Plot 7
#     fig = plt.figure()
#     plt.xlim(-16, -24)
#     plt.ylim(1, 8)
#     plt.scatter(iabs_gal[nope], nuv_i_gal[nope], color='cyan', label='no outflow pattern', s=3, alpha=1.0)
#     plt.xlabel('i-band')
#     plt.ylabel('NUV-i')
#     plt.title('I-band Absolute Magnitude vs NUV-i of Galaxies without Outflow Patterns')
#     plt.legend(loc='upper left', frameon=False, fontsize='x-small')
#
#     pdf.savefig()
#     plt.close
#
#     # Plot 7
#     fig = plt.figure()
#     plt.xlim(-16, -24)
#     plt.ylim(1, 8)
#     plt.scatter(iabs_gal, nuv_i_gal, color='lightsteelblue', label='all galaxies', s=3, alpha=.8)
#     plt.scatter(iabs_gal[nope], nuv_i_gal[nope], color='cyan', label='no outflow pattern', s=3, alpha=1.0)
#     plt.scatter(iabs_gal[clear], nuv_i_gal[clear], color='red', marker='^', label='Clear outflow', s=6 , alpha=1)
#     # plt.scatter(iabs_gal[unclear], nuv_i_gal[unclear], color='lightsteelblue', label='unclear outflow', s=5, alpha=1.0)
#     plt.xlabel('i-band')
#     plt.ylabel('NUV-i')
#     plt.title('I-band Absolute Magnitude vs NUV-i')
#     plt.legend(loc='upper left', frameon=False, fontsize='x-small')
#
#     pdf.savefig()
#     plt.close
# #
# #
#
os.system("open %s &" % filename)
