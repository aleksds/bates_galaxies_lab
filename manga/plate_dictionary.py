# Rebecca Minsley
# Tuesday July 10th, 2018

import os
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
from astropy import units as u
from os.path import isdir, join
import csv
import pandas as pd
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


mpl8_dir = os.environ['MANGADIR_MPL8']  # directory should be saved in .bash_profile

drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')    # read in information from drpall file
drpdata = drp[1].data

nsa = fits.open(mpl8_dir + '1-nsa_v1_0_1.fits')     # read in information from NSA catalog
nsa_data = nsa[1].data

blah = drpdata.field('srvymode') == 'MaNGA dither'  # boolian where True = Manga Galaxies

stellar_mass = drpdata.field('nsa_sersic_mass')
ba = drpdata.field('nsa_sersic_ba')[blah]           # sersic profile of Manga Galaxies
plateifu = drpdata.field('plateifu')

manga_ra = drpdata.field('objra')
manga_dec = drpdata.field('objdec')
nsa_ra = nsa_data.field('RA')
nsa_dec = nsa_data.field('DEC')

# read in MPA-JHU catalog information
gal = fits.open(mpl8_dir + 'gal_info_dr7_v5_2.fit')
galdata = gal[1].data
sfr = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')
sfrdata = sfr[1].data
mass = fits.open(mpl8_dir + 'gal_totsfr_dr7_v5_2.fits')
massdata = mass[1].data


# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# match the two catalogs on the basis of right ascension and declination
c = SkyCoord(ra=manga_ra*u.degree, dec=manga_dec*u.degree)
catalog = SkyCoord(ra=nsa_ra*u.degree, dec=nsa_dec*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)

# match the MPA-JHU and MaNGA catalogs on the basis of RA and DEC
mpa_ra = galdata.field('RA')
mpa_dec = galdata.field('DEC')
bad = mpa_dec < -90
mpa_dec[bad] = -90
mpa_cat = SkyCoord(ra=mpa_ra*u.degree, dec=mpa_dec*u.degree)
idx_mpa, d2d_mpa, d3d_mpa = c.match_to_catalog_sky(mpa_cat)


# define th50, th90, and concentration arrays
th50_manga = nsa_data.field('ELPETRO_TH50_R')[idx]
th90_manga = nsa_data.field('ELPETRO_TH90_R')[idx]
c_manga = (th90_manga / th50_manga)[blah]

sfr_manga = sfrdata.field('MEDIAN')[idx_mpa]
mass_manga = massdata.field('MEDIAN')[idx_mpa]


# define edge-on, late-type galaxies
late = c_manga < 2.6
ba_gal_late = ba[late]
edge = abs(ba_gal_late) < 0.3
ba_gal_late_edge = ba_gal_late[edge]
sorted_index = np.argsort(plateifu[blah][late][edge])

# mass and starformation rate
sfr_manga_gal = sfr_manga[blah][late][edge]
sfr_manga_gal = sfr_manga_gal[sorted_index]
mass_manga_gal = mass_manga[blah][late][edge]
mass_manga_gal = mass_manga_gal[sorted_index]

good_plates = np.sort(plateifu[blah][late][edge])
total = len(good_plates)
print("total late-type edge on galaxies " + str(total))


# empty arrays for data frame
re = np.zeros(total)
vel_rotation = np.zeros(total)
xig = np.zeros(total)
xi = np.zeros(total)
eta = np.zeros(total)

# creating dicitonary
for i in range(total):
    plate, ifu = good_plates[i].split('-')
    name = mpl8_dir + 'HYB10-MILESHC-MILESHC/' + str(plate) + '/' + \
               str(ifu) + '/manga-' + str(plate) + '-' + \
               str(ifu) + '-MAPS-HYB10-MILESHC-MILESHC.fits.gz'

    if os.path.isfile(name):
        match = np.where((drpdata.field('plate') == int(plate)) & (drpdata.field('ifudsgn') == str(ifu)))
        hdu_dap = fits.open(name)
        print(str(plate) + '-' + str(ifu))

        # pulling data for specified galaxy
        emlc = channel_dictionary(hdu_dap, 'EMLINE_SFLUX')
        dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Ha-6564'], :, :]
        dap_ha_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[emlc['Ha-6564'], :, :]

        # read in the position angle and axis from from the NSA
        theta = drpdata.field('nsa_sersic_phi')[match][0]
        ba = drpdata.field('nsa_sersic_ba')[match][0]
        inc = np.arccos(ba)

        # creating map
        # create arrays that correspond to x and y coordinates for each spaxel
        size = len(hdu_dap['EMLINE_GFLUX'].data[0, :])
        xpos = np.zeros([size, size])
        ypos = np.zeros([size, size])
        for j in range(0, size):
            for k in range(0, size):
                xpos[j, k] = k
                ypos[j, k] = j

        # redefine x=0 and y=0 to be in the center of the field
        xprime = xpos - np.median(xpos)
        yprime = ypos - np.median(ypos)

        # compute the on-sky x and y coordiates defined by the major axis
        trad = np.radians(theta - 90)
        xproj = xprime * np.cos(trad) + yprime * np.sin(trad)
        yproj = xprime * np.sin(trad) * (-1.) + yprime * np.cos(trad)
        zproj = yproj / np.sin(inc) * np.cos(inc)

        # calculate the radius of each pixel in the plane of the disk [units: pixels]
        radpix = np.sqrt(xproj ** 2 + (yproj / ba) ** 2)

        # figure out the conversion between pixel size and kpc
        z = drpdata.field('nsa_z')[match][0]
        radkpc = radpix / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

        # compute values along the x and y axis of the disk and the z axis above the disk
        xproj_kpc = xproj / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
        yproj_kpc = yproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
        zproj_kpc = zproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
        cen = abs(xproj_kpc) < (1.1 * u.kpc)

        radkpc_map = np.zeros([size, size])
        xproj_kpc_map = np.zeros([size, size])
        yproj_kpc_map = np.zeros([size, size])
        zproj_kpc_map = np.zeros([size, size])

        for j in range(0, size):
            for k in range(0, size):
                if dap_ha_sivar[j, k] > 0.:
                    radkpc_map[j, k] = radkpc[j, k] / u.kpc
                    yproj_kpc_map[j, k] = yproj_kpc[j, k] / u.kpc
                    zproj_kpc_map[j, k] = zproj_kpc[j, k] / u.kpc
                else:
                    radkpc_map[j, k] = radkpc[j, k] / u.kpc * (-1.)

        for j in range(0, 1):
            # DAP: emission-line quantities
            dap_vel = hdu_dap['EMLINE_GVEL'].data[j, :, :]
            dap_vel_ivar = hdu_dap['EMLINE_GVEL_IVAR'].data[j, :, :]
            dap_vel_err = np.sqrt(1. / dap_vel_ivar)
            dap_sig = hdu_dap['EMLINE_GSIGMA'].data[j, :, :]
            dap_sig_ivar = hdu_dap['EMLINE_GSIGMA'].data[j, :, :]
            dap_sig_err = np.sqrt(1. / dap_sig_ivar)
            dap_flux = hdu_dap['EMLINE_GFLUX'].data[j, :, :]
            dap_mask = hdu_dap['EMLINE_GFLUX_MASK'].data[j, :, :]
            dap_ivar = hdu_dap['EMLINE_GFLUX_IVAR'].data[j, :, :]
            dap_err = np.sqrt(1. / dap_ivar)
            dap_sflux = hdu_dap['EMLINE_SFLUX'].data[j, :, :]
            dap_smask = hdu_dap['EMLINE_SFLUX_MASK'].data[j, :, :]
            dap_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[j, :, :]

            # H_alpha/H_beta
            dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Ha-6564'], :, :]
            dap_hb_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Hb-4862'], :, :]
            ha_hb = dap_ha_sflux / dap_hb_sflux

            # calculating noise variables
            dap_serr = np.sqrt(1. / dap_sivar)
            s_n_ha = dap_ha_sflux * np.sqrt(dap_ha_sivar)    # signal to noise h alpha

            # creating maps flipped across the major axis
            vel_flip = np.zeros([size, size])
            ivar_flip = np.zeros([size, size])
            snha_flip = np.zeros([size, size])
            for m in range(0, size):
                for n in range(0, size):
                    dist = np.sqrt((xproj_kpc - xproj_kpc[m, n]) ** 2 + (yproj_kpc + yproj_kpc[m, n]) ** 2)
                    test = (dist == np.min(dist))
                    vel_flip[m, n] = dap_vel[test]
                    ivar_flip[m, n] = dap_vel_ivar[test]
                    snha_flip[m, n] = s_n_ha[test]

            # identifying bad spaxels
            snb = s_n_ha < 5
            snb_flip = snha_flip < 5

            # calculating asymmetry parameter including spaxels with s_n(h alpha) < 5
            delta_vel = dap_vel-vel_flip
            delta_vel[snb] = 0
            delta_vel[snb_flip] = 0
            denom = np.sqrt((1/dap_ivar)+(1/ivar_flip))

            re_arcsec = drpdata.field('nsa_elpetro_th50_r')[match] + np.array([1])
            re_kpc = re_arcsec * u.arcsec / cosmo.arcsec_per_kpc_proper(z)
            cen = abs(radkpc_map) < (re_kpc / u.kpc)
            re_large = np.ravel(radkpc) > re_kpc
            re_large_positive = np.ravel(yproj_kpc_map)[re_large] > 0
            re_large_negative = np.ravel(yproj_kpc_map)[re_large] < 0

            # calculating asymmetry parameter removing the spaxels with s_n(h alpha) < 5
            frac_g = delta_vel/denom
            ha_good_large_positive = np.ravel(s_n_ha)[re_large][re_large_positive] > 5
            ha_good_large_negative = np.ravel(s_n_ha)[re_large][re_large_negative] > 5
            poop_positive_good = (np.std(np.ravel(frac_g)[re_large][re_large_positive][ha_good_large_positive]))
            poop_negative_good = (np.std(np.ravel(frac_g)[re_large][re_large_negative][ha_good_large_negative]))
            poop_tmp_good = (poop_positive_good + poop_negative_good) / 2.      # asymmetry parameter value with bad spaxels removed
            print(poop_tmp_good)

            # calculating eta or the velocity dispersion to rotation variable
            # velocity dispersion
            ha_good_large = np.ravel(s_n_ha)[re_large] > 5
            med_vel_disp = np.median(np.ravel(dap_sig)[re_large][ha_good_large])
            print('the velocity dispersion is', med_vel_disp)

            # velocity roation with ave of vrot
            major = abs(yproj_kpc_map) < 1.0  # within 1 kpc of major axis
            vel_major = (dap_vel)[major]  # might need np.ravel() ?

            vrot = np.percentile(vel_major, 90)
            vrot_ave = np.round(np.mean(vrot, dtype=None), 5)
            print('the rotational velocity is', vrot_ave)
            # Vrot/Vdisp
            med_vel_disp_vrot_ave = med_vel_disp / vrot_ave
            print('the velocity dipsersion to rotation ratio is', med_vel_disp_vrot_ave)

            # arrays for data frame
            re[i] = re_kpc.value[0]
            vel_rotation[i] = vrot
            xig[i] = poop_tmp_good
            eta[i] = med_vel_disp_vrot_ave


# calculating starformation density
sf_density = sfr_manga_gal / ((re ** 2) * math.pi)

# calculating mass from vel rotation v^2 = GM/r => M = (v^2 *r)/G in SOLAR MASS
vel_mass = ((((vel_rotation * (u.km / u.s)) ** 2) * re_kpc/G).to('Msun')).value

data = good_plates, re, vel_rotation, mass_manga_gal, vel_mass, sfr_manga_gal, sf_density, xig, eta
df = pd.DataFrame({'Galaxy': data[0], 're (kpc)': data[1], 'rotational velocity': data[2],
                   'rotational velocity mass (solar mass)': data[3], 'star formation rate': data[4],
                   'star formation density': data[5], 'assymetry': data[6], 'eta average vrot': data[7]})
df.to_csv('galaxy_dictionary.csv')

filename = 'mpl8_sfr_density'
with PdfPages(filename) as pdf:
    plt.scatter(mass_manga_gal, sfr_manga_gal, s=0.1, color='lightsteelblue', marker='s')
    # plt.xlim(8, 12)
    # plt.ylim(-2, 1.5)
    plt.xlabel('Log(Stellar Mass)')
    plt.ylabel('Log(Star Formation Rate)')
    pdf.savefig()
    plt.close()
os.system("open %s &" % filename)



