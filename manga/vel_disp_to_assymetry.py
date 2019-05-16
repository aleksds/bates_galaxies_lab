# Rebecca Minsley
# Tuesday July 10th, 2018

import os
import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from os.path import isdir, join
import glob
import csv
import math
import pickle


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

# minimum and maximum emission-line fluxes for plot ranges
fmin = 1e-19
fmax = 1e-16

def daplot(quantity, qmin, qmax):
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    plt.imshow(quantity, origin='lower', interpolation='nearest',
               norm=colors.LogNorm(vmin=qmin, vmax=qmax), cmap=cm.coolwarm)




# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

mpl7_dir = os.environ['MANGADIR_MPL7']  # Be aware that this directory syntax might need revision

# read in information from drpall file
drp = fits.open(mpl7_dir + 'drpall-v2_4_3.fits')
drpdata = drp[1].data

# read in information from NSA catalog
nsa = fits.open(mpl7_dir + '1-nsa_v1_0_1.fits')
nsa_data = nsa[1].data

plate = 7977
galaxy = 12704
match = np.where((drpdata.field('plate') == plate) & (drpdata.field('ifudsgn') == str(galaxy)))
mangaid = drpdata[match].field('mangaid')[0]
nsaid = drpdata[match].field('nsa_nsaid')[0]
yo = np.where(nsa_data.field('NSAID') == nsaid)
name = nsa_data[yo[0][0]].field('IAUNAME')
what = drpdata[match].field('nsa_iauname')[0]

ba = drpdata.field('nsa_sersic_ba')
plateifu = drpdata.field('plateifu')
manga_ra = drpdata.field('objra')
manga_dec = drpdata.field('objdec')
nsa_ra = nsa_data.field('RA')
nsa_dec = nsa_data.field('DEC')

# match the two catalogs on the basis of right ascension and declination
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
# update the concentration and ba arrays
c_manga_gal = c_manga[blah]
ba_gal = ba[blah]
# define edge-on, late-type galaxies
late = c_manga_gal < 2.6
ba_gal_late = ba_gal[late]
edge = abs(ba_gal_late) < 0.3
ba_gal_late_edge = ba_gal_late[edge]
good_plates = np.sort(plateifu[blah][late][edge])
elp_50 = th50_manga[blah][late][edge]

total = len(good_plates)
print("total late-type edge on galaxies " + str(total))
xib = np.zeros(total)
xig = np.zeros(total)
poo_n = np.zeros(total)


filename = 'xi_plotz.pdf'
with PdfPages(filename) as pdf:
    for i in range(total):
        plate, ifu = good_plates[i].split('-')
        name = mpl7_dir + 'HYB10-GAU-MILESHC/' + str(plate) + '/' + \
               str(ifu) + '/manga-' + str(plate) + '-' + \
               str(ifu) + '-MAPS-HYB10-GAU-MILESHC.fits.gz'

        if os.path.isfile(name):
            match = np.where((drpdata.field('plate') == int(plate)) & (drpdata.field('ifudsgn') == str(ifu)))
            hdu_dap = fits.open(name)
            print(str(plate) + '-' + str(ifu))

            # pulling data for specified galaxy
            # hdu_dap.info()
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
                    if (dap_ha_sivar[j, k] > 0.):
                        radkpc_map[j, k] = radkpc[j, k] / u.kpc
                        yproj_kpc_map[j, k] = yproj_kpc[j, k] / u.kpc
                        zproj_kpc_map[j, k] = zproj_kpc[j, k] / u.kpc
                    else:
                        radkpc_map[j, k] = radkpc[j, k] / u.kpc * (-1.)

            # hdu_dap['EMLINE_GFLUX'].header
            # C01     = 'OII-3727'           / Data in channel 1
            # C02     = 'OII-3729'           / Data in channel 2
            # C03     = 'Hthe-3798'          / Data in channel 3
            # C04     = 'Heta-3836'          / Data in channel 4
            # C05     = 'NeIII-3869'         / Data in channel 5
            # C06     = 'Hzet-3890'          / Data in channel 6
            # C07     = 'NeIII-3968'         / Data in channel 7
            # C08     = 'Heps-3971'          / Data in channel 8
            # C09     = 'Hdel-4102'          / Data in channel 9
            # C10     = 'Hgam-4341'          / Data in channel 10
            # C11     = 'HeII-4687'          / Data in channel 11
            # C12     = 'Hb-4862 '           / Data in channel 12
            # C13     = 'OIII-4960'          / Data in channel 13
            # C14     = 'OIII-5008'          / Data in channel 14
            # C15     = 'HeI-5877'           / Data in channel 15
            # C16     = 'OI-6302 '           / Data in channel 16
            # C17     = 'OI-6365 '           / Data in channel 17
            # C18     = 'NII-6549'           / Data in channel 18
            # C19     = 'Ha-6564 '           / Data in channel 19
            # C20     = 'NII-6585'           / Data in channel 20
            # C21     = 'SII-6718'           / Data in channel 21
            # C22     = 'SII-6732'           / Data in channel 22

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
                s_n_ha = dap_ha_sflux * np.sqrt(dap_ha_sivar)     # signal to noise h alpha


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
                poop_tmp_good = (poop_positive_good + poop_negative_good) / 2.  # asymmetry parameter value with bad spaxels removed
                xig[i] = poop_tmp_good
                print(poop_tmp_good)

                # calculating eta or the velocity dispersion to rotation variable
                # velocity dispersion
                ha_good_large = np.ravel(s_n_ha)[re_large] > 5
                med_vel_disp = np.median(np.ravel(dap_sig)[re_large][ha_good_large])
                print('the velocity dispersion is', med_vel_disp)

                # velocity roation
                major = yproj_kpc_map < 1.0  # within 1 kpc of major axis
                vel_major = (dap_vel)[major]  # might need np.ravel() ?
                vrot = np.percentile(vel_major, 90)  # 90th percentile make sense?
                vrot_max = np.round(np.mean(vrot, dtype=None), 5)
                print('the rotational velocity is', vrot_max)
                # Vrot/Vdisp
                med_vel_disp_vrot_max = med_vel_disp / vrot_max
                print('the velocity dipsersion to rotation ratio is', med_vel_disp_vrot_max)

                if xig[i] > 2 and med_vel_disp_vrot_max > 0.3:
                    print(str(ha_good_large_positive))
                    # plot 1: ha flux
                    fig = plt.figure()
                    # plt.title(str(good_plates(i)) + '\n' + 'η50' + str(med_vel_disp_vrot_max))

                    plt()
                    ax = fig.add_subplot(2, 4, 1)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow((dap_ha_sflux *1.e-17), origin='lower',interpolation='nearest',
                               norm=colors.LogNorm(vmin=fmin, vmax=fmax), cmap=cm.coolwarm)
                    plt.title(str(good_plates[i]) + r'$ H\alpha$', fontsize=8)
                    cb = plt.colorbar(fraction=0.047, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)

                    # # plot 1: ha flux
                    # dap_ha_sflux[cen] = 0
                    # ax = fig.add_subplot(2, 4, 5)
                    # ax.patch.set_facecolor('white')
                    # ax.set_xlim(0, size)
                    # ax.set_ylim(0, size)
                    # ax.set_xticklabels(())
                    # ax.set_yticklabels(())
                    # plt.imshow((dap_ha_sflux *1.e-17), origin='lower',interpolation='nearest',
                    #            norm=colors.LogNorm(vmin=fmin, vmax=fmax), cmap=cm.coolwarm)
                    # plt.title(r'$H\alpha$', fontsize=8)
                    # # cb = plt.colorbar(shrink=.20, pad=.03, orientation='horizontal')
                    # # cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=5)


                    # plot 2: Ha/Hb emission
                    ax = fig.add_subplot(2, 4, 2)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(ha_hb, origin='lower', interpolation='nearest',
                               norm=colors.LogNorm(vmin=0.1, vmax=10), cmap=cm.coolwarm)
                    plt.title(r'$H\alpha / H\beta$', fontsize=8)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)

                    # plot 6: masked Ha/Hb emission
                    ha_hb[cen] = 0
                    ax = fig.add_subplot(2, 4, 6)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(ha_hb, origin='lower', interpolation='nearest',
                               norm=colors.LogNorm(vmin=0.1, vmax=10), cmap=cm.coolwarm)
                    plt.title(r'$H\alpha / H\beta$', fontsize=8)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)


                    # plot 3: velocity
                    ax = fig.add_subplot(2, 4, 3)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_vel, origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm, vmin=-250, vmax=250)
                    plt.title('Ha Vel', fontsize=8)
                    plt.title(r'$v_{gas}$', verticalalignment='center', fontsize=8)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)


                    # plot 7: masked velocity
                    dap_vel[cen] = 0
                    ax = fig.add_subplot(2, 4, 7)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_vel, origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm, vmin=-250, vmax=250)
                    # plt.colorbar()
                    plt.title('Ha Vel', fontsize=8)
                    plt.title(r'$v_{gas}$', verticalalignment='center', fontsize=8)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)


                    # plot 4: assymetry
                    ax = fig.add_subplot(2, 4, 4)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_yticklabels(())
                    ax.set_xticklabels(())
                    plt.tick_params(axis='both',length=0)
                    plt.imshow(frac_g, origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm, vmin=-40, vmax=40)
                    plt.title('asymmetry', verticalalignment='center', fontsize=8)
                    ax.tick_params(axis='both', reset=False, which='both', length=0, width=0)
                    plt.text(1, -8., r'$ξ_{sn>5} =$' + ' ' + str(round(poop_tmp_good,4)), fontsize=7)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)


                    # plot 8: v - v_flip
                    delta_vel[cen] = 0
                    ax = fig.add_subplot(2, 4, 8)
                    ax.patch.set_facecolor('white')
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_yticklabels(())
                    ax.set_xticklabels(())
                    plt.tick_params(axis='both',length=0)
                    plt.imshow(delta_vel, origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm, vmin=-40, vmax=40)
                    plt.title(r'$v - v_{flip}$', verticalalignment='center', fontsize=8)
                    ax.tick_params(axis='both', reset=False, which='both', length=0, width=0)
                    plt.text(1, -8., r'$ξ_{sn>5} =$' + ' ' + str(round(poop_tmp_good,4)), fontsize=7)
                    cb = plt.colorbar(fraction=0.046, pad=0.02)
                    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=3)







                    pdf.savefig(facecolor=fig.get_facecolor(), edgecolor='none')
                    plt.close()


# with open('kinematic_var.pickle', 'wb') as handle:
#     pickle.dump((xib,xig), handle, protocol=pickle.HIGHEST_PROTOCOL)
os.system("open %s &" % filename)








