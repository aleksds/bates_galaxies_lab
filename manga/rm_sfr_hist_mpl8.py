import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv

# MPL 8 data
mpl8_dir = os.environ['MANGADIR_MPL8']  # Be aware that this directory syntax might need revision
drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')    # read in information from drpall file
drpdata = drp[1].data
ba = drpdata.field('nsa_sersic_ba')
stellar_mass = drpdata.field('nsa_sersic_mass')
absmag = drpdata.field('nsa_elpetro_absmag')
gr = absmag[:, 3] - absmag[:, 4]
filename = 'histograms.pdf'
plateifu = drpdata.field('plateifu')

# read in MPA-JHU catalog information
gal = fits.open(mpl8_dir + 'gal_info_dr7_v5_2.fit')
galdata = gal[1].data
sfr = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')
sfrdata = sfr[1].data
mass = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')
massdata = mass[1].data

# read in information from NSA catalog
nsa = fits.open(mpl8_dir + '1-nsa_v1_0_1.fits')     # read in information from NSA catalog
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
what = drpdata[match].field('nsa_iauname')[0]
print(what)

# match the NSA and MaNGA catalogs on the basis of right ascension and declination
manga_ra = drpdata.field('objra')
manga_dec = drpdata.field('objdec')
nsa_ra = nsa_data.field('RA')
nsa_dec = nsa_data.field('DEC')
c = SkyCoord(ra=manga_ra * u.degree, dec=manga_dec * u.degree)
catalog = SkyCoord(ra=nsa_ra * u.degree, dec=nsa_dec * u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)

# match the MPA-JHU and MaNGA catalogs on the basis of RA and DEC
mpa_ra = galdata.field('RA')
mpa_dec = galdata.field('DEC')
bad = mpa_dec < -90
mpa_dec[bad] = -90
mpa_cat = SkyCoord(ra=mpa_ra * u.degree, dec=mpa_dec * u.degree)
idx_mpa, d2d_mpa, d3d_mpa = c.match_to_catalog_sky(mpa_cat)

# define th50, th90, and concentration arrays
th50 = nsa_data.field('ELPETRO_TH50_R')
th90 = nsa_data.field('ELPETRO_TH90_R')

th50_manga = th50[idx]
th90_manga = th90[idx]
c_manga = th90_manga / th50_manga

sfr_manga = sfrdata.field('MEDIAN')[idx_mpa]
mass_manga = massdata.field('MEDIAN')[idx_mpa]

# get just the MaNGA galaxies
blah = drpdata.field('srvymode') == 'MaNGA dither'
print(len(c_manga[blah]))

# update the concentration and ba arrays
c_manga_gal = c_manga[blah]
ba_gal = ba[blah]
stellar_mass_gal = stellar_mass[blah]
gr_gal = gr[blah]
sfr_manga_gal = sfr_manga[blah]
mass_manga_gal = mass_manga[blah]

# find late-type manga galaxies
late = c_manga_gal < 2.6

# update concentration and ba arrays
ba_gal_late = ba_gal[late]
c_manga_gal_late = c_manga_gal[late]
stellar_mass_gal_late = stellar_mass_gal[late]
gr_gal_late = gr_gal[late]
sfr_manga_gal_late = sfr_manga_gal[late]
mass_manga_gal_late = mass_manga_gal[late]

# find edge-one late-type galaxies
edge = abs(ba_gal_late) < 0.3

# update concentration and ba arrays
c_manga_gal_late_edge = c_manga_gal_late[edge]
ba_gal_late_edge = ba_gal_late[edge]
stellar_mass_gal_late_edge = stellar_mass_gal_late[edge]
gr_gal_late_edge = gr_gal_late[edge]
sfr_manga_gal_late_edge = sfr_manga_gal_late[edge]
mass_manga_gal_late_edge = mass_manga_gal_late[edge]



print(len(c_manga_gal_late_edge))
print(len(ba_gal_late_edge))
print(np.sort(plateifu[blah][late][edge]))
print(len(plateifu[blah][late][edge]))
###### Importing data from csv file
# Empty lists for csv file data
outflow_assign = []                                                     # for the data on the file about outflow
br_assign = []                                                          # for the data on the file about b/r sep

# open file containing classification of late-type-edge-on galaxies
# with open('em_outflow.csv', 'r') as csv_file:
#     for row in csv.reader(csv_file, delimiter=','):
#         outflow_assign.append(row[3])
#         br_assign.append(row[4])
#
# del outflow_assign[0]                                                    # deleting the column label outflow the csv file data set
# del br_assign[0]                                                     # deleting the column label blue/red sep from the csv file data set

# Empty lists for categorization
c1 = np.zeros(445)                                                           # clear outflow
c2 = np.zeros(445)                                                           # not clear outflow

# counters to see how many galaxies lie in each outflow category
g = 0                                                                  # category 1 counter
n = 0                                                                  # category 2 counter
# outflows category assignment
# for i in range(0,445):
#     if outflow_assign[i] == 'yes':
#         c1[i] = 10
#         g = g + 1
#     else:
#         c2[i] = -10
#         n = n + 1
# print('outflow # = '+ str(g))

# outflow categories bool arrays
outflow = np.array(c1) > 0
no_flow = np.array(c2) < 0

#ouflow arrays
stellar_flow = stellar_mass_gal_late_edge[outflow]



with PdfPages(filename) as pdf:
    # plot 1: b/a ratio
    fig = plt.figure()
    # plt.xlim(xmin=0.0, xmax=1.05)
    # plt.hist(ba_gal, color='lightsteelblue', bins=np.arange(0.05, 1.02, .02), label='All MPL-7 galaxies', alpha=.4)
    # plt.hist(ba_gal_late_edge, color='navy', bins=np.arange(0.05, 1.02, .02), label='Late-type-edge-on Galaxies', alpha=.4)
    # plt.hist(ba_gal_late_edge[outflow], color='red', bins=np.arange(0.05, 1.02, .02), label='Outflows', alpha=.7)
    # plt.xlabel('b/a ratio')
    # plt.ylabel('Number of Galaxies')
    # plt.title('Distribution of b/a Values in MaNGA sample')
    # plt.legend(loc='upper left', frameon=False, fontsize='x-small')
    #
    # pdf.savefig()
    # plt.close

    log_stellar_mass_gal = np.log10(stellar_mass_gal)
    use = np.isfinite(log_stellar_mass_gal)

    log_stellar_mass_gal_late_edge = np.log10(stellar_mass_gal_late_edge)
    yo = np.isfinite(log_stellar_mass_gal_late_edge)

    log_stellar_mass_gal_late = np.log10(stellar_mass_gal_late)
    mo = np.isfinite(log_stellar_mass_gal_late)

    log_stellar_flow = np.log10(stellar_flow)
    flo = np.isfinite(log_stellar_flow)

    # plot 2: stellar mass
    fig = plt.figure()
    plt.hist(log_stellar_mass_gal[use], color='lightsteelblue', bins=np.arange(8.0, 12., .1), label='All MPL-7 galaxies', alpha=.4)
    plt.hist(log_stellar_mass_gal_late[mo], color='cyan', bins=np.arange(8.0, 12., .1), label='Late-type Galaxies', alpha=.4)
    plt.hist(log_stellar_mass_gal_late_edge[yo], color='navy', bins=np.arange(8.0, 12., .1),label='Late-type-edge-on Galaxies', alpha=.4)
    plt.hist(log_stellar_flow[flo], color='red', bins=np.arange(8.0, 12., .1),label='Outflows', alpha=.7)
    plt.xlim(8, 12)
    plt.xlabel('Stellar Mass (logM☉)')
    plt.ylabel('Number of Galaxies')
    plt.title('Distribution of Stellar Mass in MaNGA sample')
    plt.legend(loc='upper left', frameon=False, fontsize='x-small')

    pdf.savefig()
    plt.close

    # # plot 3: g-r
    # fig = plt.figure()
    # plt.hist(gr_gal, bins=np.arange(0.0, 1.0, 0.05), color='lightsteelblue', label='All MPL-7 galaxies', alpha=.4)
    # plt.hist(gr_gal_late, bins=np.arange(0.0, 1.0, 0.05), color='cyan', label='Late-type Galaxies', alpha=.4)
    # plt.hist(gr_gal_late_edge, bins=np.arange(0.0, 1.0, 0.05), color='navy', label='Late-type-edge-on Galaxies', alpha=.4)
    # plt.hist(gr_gal_late_edge[outflow], bins=np.arange(0.0, 1.0, 0.05), color='red', label='Outflow', alpha=.7)
    # plt.xlabel('g-r')
    # plt.ylabel('Number of Galaxies')
    # plt.title('Distribution of g-r in MaNGA sample')
    # plt.legend(loc='upper left', frameon=False, fontsize='x-small')
    #
    # pdf.savefig()
    # plt.close

    # # concentration
    # fig = plt.figure()
    # plt.xlim(xmin=0, xmax=4.5)
    # plt.hist(c_manga_gal, bins=np.arange(.25, 4.1, 0.1), color='lightsteelblue', label='All MPL-7 galaxies', alpha=.4)
    # plt.hist(c_manga_gal_late_edge, bins=np.arange(.25, 4.1, 0.1), color='navy', label='Late-type-edge-on Galaxies', alpha=.4)
    # plt.hist(c_manga_gal_late_edge[outflow], bins=np.arange(.25, 4.1, 0.1), color='red', label='Outflow', alpha=.7)
    # plt.xlabel('Concentration')
    # plt.ylabel('Number of Galaxies')
    # plt.title('Distribution of Concentration in MaNGA sample')
    # plt.legend(loc='upper left', frameon=False, fontsize='x-small')
    #
    # pdf.savefig()
    # plt.close
    #
    # fig = plt.figure()
    #
    plt.scatter(mass_manga_gal, sfr_manga_gal, s=0.1, color='lightsteelblue', marker='s')
    #plt.scatter(mass_manga_gal_late, sfr_manga_gal_late, s=0.1, color='blue', marker='o')
    plt.scatter(mass_manga_gal_late_edge, sfr_manga_gal_late_edge, s=0.5, marker='*', color='navy')
    # plt.scatter(mass_manga_gal_late_edge[outflow], sfr_manga_gal_late_edge[outflow], s=0.5, marker='*', color='red')
    plt.xlim(8, 12)
    plt.ylim(-2, 1.5)
    plt.xlabel('Log(Stellar Mass)')
    plt.ylabel('Log(Star Formation Rate)')
    pdf.savefig()
    plt.close()

os.system("open %s &" % filename)

