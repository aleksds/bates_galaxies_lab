import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.colors as colors
from matplotlib import cm
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from os.path import isdir, join
import glob
from marvin.utils.general.images import showImage



# function for plotting maps of relevant quantities
def daplot(quantity, qmin, qmax):
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    plt.imshow(quantity, origin='lower', interpolation='nearest',
               norm=colors.LogNorm(vmin=qmin, vmax=qmax), cmap=cm.coolwarm)
    plt.colorbar()


def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:]) - 1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


# minimum and maximum emission-line fluxes for plot ranges
fmin = 1e-19
fmax = 1e-16

# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

mpl7_dir = os.environ['MANGADIR_MPL7']  # Be aware that this directory syntax might need revision
drp = fits.open(mpl7_dir + 'drpall-v2_4_3.fits')
drpdata = drp[1].data

ba = drpdata.field('nsa_sersic_ba')
plateifu = drpdata.field('plateifu')

# read in information from NSA catalog
nsa = fits.open(mpl7_dir + '1-nsa_v1_0_1.fits')
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

# update the concentration and ba arrays
c_manga_gal = c_manga[blah]
ba_gal = ba[blah]

# define edge-on, late-type galaxies
late = c_manga_gal < 2.6
ba_gal_late = ba_gal[late]
edge = abs(ba_gal_late) < 0.3
ba_gal_late_edge = ba_gal_late[edge]
print('there are', len(ba_gal_late_edge),'edge-on and late-type galxies in the mpl7 survey')
good_plates = np.sort(plateifu[blah][late][edge])


filename = 'all_good_plots.pdf'
with PdfPages(filename) as pdf:
    for i in range(0, 3):  # len(good_plates)):
        plate, ifu = good_plates[i].split('-')
        name = mpl7_dir + 'HYB10-GAU-MILESHC/' + str(plate) + '/' + str(ifu) + '/manga-' + str(plate) + '-' + str(
            ifu) + '-MAPS-HYB10-GAU-MILESHC.fits.gz'

        if os.path.isfile(name):
            match = np.where((drpdata.field('plate') == int(plate)) & (drpdata.field('ifudsgn') == str(ifu)))
            hdu_dap = fits.open(name)

            # hdu_dap.info()
            emlc = channel_dictionary(hdu_dap, 'EMLINE_SFLUX')
            dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Ha-6564'], :, :]
            dap_hb_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Hb-4862'], :, :]
            dap_ha_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[emlc['Ha-6564'], :, :]
            dap_ha_vel = hdu_dap['EMLINE_GVEL'].data[emlc['Ha-6564'], :, :]
            dap_ha_sig = hdu_dap['EMLINE_GSIGMA'].data[emlc['Ha-6564'], :, :]


            # create arrays that correspond to x and y coordinates for each spaxel
            size = len(hdu_dap['EMLINE_GFLUX'].data[0, :])
            xpos = np.zeros([size, size])
            ypos = np.zeros([size, size])

            for j in range(0, size):
                for k in range(0, size):
                    xpos[j, k] = k
                    ypos[j, k] = j

            # read in the position angle and axis from from the NSA
            theta = drpdata.field('nsa_sersic_phi')[match][0]
            ba = drpdata.field('nsa_sersic_ba')[match][0]
            inc = np.arccos(ba)

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
            cen = abs(xproj_kpc) < (1. * u.kpc)

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

            # for j in range(0, len(hdu_dap['EMLINE_GVEL'].data)):
            # focus on the [O II] emisison line
            for j in range(0, 1):
                # print(eml[j])

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
                dap_serr = np.sqrt(1. / dap_sivar)  # *1.e-4
                dap_ha_serr = np.sqrt(1. / dap_ha_sivar)


                vel_flip = np.zeros([size, size])
                ha_hb = dap_ha_sflux / dap_hb_sflux
                ha_hb_flip =np.zeros([size,size])
                ha_vel_flip = np.zeros([size, size])

                s_n_ha = dap_ha_sflux/dap_ha_serr
                bad = (s_n_ha < 2.)
                s_n_ha[bad]=0

                s_n = dap_sflux/dap_serr
                bad = (s_n < 2.)
                s_n[bad]=0

                dap_ha_vel[bad]= -2000
                for i in range(0,len(dap_vel[i])):
                    if dap_ha_vel[i].all() < -1999:
                            dap_vel.remove(dap_vel[i].all())

                dap_vel[bad]= -2000
                for i in range(0,len(dap_vel[i])):
                    if dap_vel[i].all() < -1999:
                            dap_vel.remove(dap_vel[i].all())


                for m in range(0, size):
                    for n in range(0, size):
                        dist = np.sqrt((xproj_kpc - xproj_kpc[m, n]) ** 2 + (yproj_kpc + yproj_kpc[m, n]) ** 2)
                        test = (dist == np.min(dist))
                        vel_flip[m, n] = dap_vel[test]
                        ha_hb_flip[m, n] = ha_hb[test]
                        ha_vel_flip[m, n]= dap_ha_vel[test]


                # page 1, plot 1: galaxy coordinates in kpc
                fig = plt.figure()
                ax = fig.add_subplot(3,3,1)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(zproj_kpc_map,
                           origin='lower',
                           interpolation='nearest',
                           cmap=cm.coolwarm)
                plt.colorbar()
                plt.title('vertical distance [kpc]', fontsize=10)
                plt.suptitle('('+str(i)+')'+' '+plate+'-'+ifu)

                # page 1, plot 2: emission-line SFLUX
                ax = fig.add_subplot(3,3,2)
                daplot(dap_sflux*1.e-17, fmin, fmax)
                plt.title(' O[II] flux', fontsize=10)

    
                # page 1, plot 3: emission-line SFLUX serror
                ax = fig.add_subplot(3,3,3)
                daplot(dap_serr*1.e-17, fmin, fmax)
                plt.title('O[II] flux Error', fontsize=10)

    
                # page 1, plot 4: emission-line SFLUX signal-to-noise ratio
                ax = fig.add_subplot(3,3,4)
                daplot(s_n, 0.1, 10.)
                plt.title('signal to noise ratio of O[II] flux', fontsize=10)

    
                # page 1, plot 5: emission-line flux / Halpha flux
                ax = fig.add_subplot(3,3,5)
                daplot(dap_sflux/dap_ha_sflux, 0.1, 10.)
                plt.title('O[II] to HA ratio', fontsize=10)

    
                # page 1, plot 6: emission-line dispersion
                ax = fig.add_subplot(3,3,6)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_sig, origin='lower',
                           interpolation='nearest',
                           cmap=cm.YlOrRd, vmin=0, vmax=250)
                plt.colorbar()
                plt.title('GSIGMA', fontsize=10)

                
                # page 1, plot 7: emission-line velocity
                ax = fig.add_subplot(3,3,7)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_vel, origin='lower',
                           interpolation='nearest',
                           cmap=cm.coolwarm, vmin=-250, vmax=250)
                plt.colorbar()
                plt.title('GVEL', fontsize=10)

                # page 1, plot 8: Vflip
                ax = fig.add_subplot(3,3,8)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(vel_flip, origin='lower',
                           interpolation='nearest',
                           cmap=cm.coolwarm, vmin=-250, vmax=250)
                plt.colorbar()
                plt.title('VFLIP', fontsize=10)

                #page 1, plot 9: GVEL-Vflip
                ax = fig.add_subplot(3,3,9)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_vel-vel_flip, origin='lower',
                           interpolation='nearest',
                           cmap=cm.coolwarm, vmin=-80, vmax=80)
                plt.colorbar()
                plt.title('GVEL-VFLIP', fontsize=10)

    
                pdf.savefig()
                plt.close()

                # page 2, plot 1: plots the image of the galaxy
                fig = plt.figure()
                ax = fig.add_subplot(3, 3, 1)
                #image = showImage(plateifu=good_plates[i], show_image=False)
                #fig.suptitle(good_plates[i]+' '+ '('+str(i)+')', fontsize = 12)
                #plt.imshow(image)
                plt.axis('off')
                print(good_plates[i] + ' ' + str(i))

                # page 2, plot 2: ha flux
                ax = fig.add_subplot(3, 3, 2)
                daplot(dap_ha_sflux * 1.e-17, fmin, fmax)
                plt.title('Ha flux', fontsize=8)

                # page 2, plot 3: ha error
                ax = fig.add_subplot(3, 3, 3)
                daplot(dap_serr * 1.e-17, fmin, fmax)
                plt.title('Ha flux Error', fontsize=8)

                # page 2, plot 4: s/n
                ax = fig.add_subplot(3, 3, 4)
                daplot(s_n_ha, 0.1, 10.)
                plt.title('signal to noise ratio of Halpha flux', fontsize=8)

                # page 2, plot 5: Ha/Hb emission
                ax = fig.add_subplot(3, 3, 5)
                daplot(dap_ha_sflux / dap_hb_sflux, 0.1, 10.)
                plt.title('H alpha / H beta ratio', fontsize=8)


                # page 2, plot 6: Ha/Hb - Ha/Hb flip
                ax = fig.add_subplot(3, 3, 6)
                daplot(ha_hb - ha_hb_flip, 0.1, 10.)
                plt.title('Ha/Hb - Ha/Hb flip', fontsize=8)

                # page 1, plot 7: Ha emission-line dispersion
                ax = fig.add_subplot(3,3,7)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_sig, origin='lower',
                           interpolation='nearest',
                           cmap=cm.YlOrRd, vmin=0, vmax=250)
                plt.colorbar()
                plt.title('Ha GSIGMA', fontsize=10)

                # page 2, plot 8: Halpha velovity
                ax = fig.add_subplot(3, 3, 8)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_ha_vel, origin='lower',
                           interpolation='nearest',
                           cmap=cm.seismic, vmin=-250, vmax=250)
                plt.colorbar()
                plt.title('Ha Velocity', fontsize=8)


                #plot 9: Halpha vel - halpha vel flipped
                ax = fig.add_subplot(3, 3, 9)
                ax.set_xlim(0, size)
                ax.set_ylim(0, size)
                ax.set_xticklabels(())
                ax.set_yticklabels(())
                plt.imshow(dap_ha_vel - ha_vel_flip, origin='lower',
                           interpolation='nearest',
                           cmap=cm.seismic, vmin=-250, vmax=250)
                plt.colorbar()
                plt.title('Ha Vel - ha Vel flipped', fontsize=8)


                pdf.savefig()
                plt.close()

os.system("open %s &" % filename)




