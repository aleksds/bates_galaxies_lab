import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
import matplotlib.colors as colors
from matplotlib import cm
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from os.path import isdir, join
import glob
import csv
import math

# Categories of Galaxies Outflows
#   1) clear outflow
#     1*) exciting outflow along major axis, inflow along minor axis
#   2) clear opposite outflow
#   3) no clear blue/red extraplanar separation
#
#   0) clear blue/red separation but no clear orientation

# EMLINE_GVEL_IVAR gives inverse variance = 1 / sigma^2 (where sigma is the 1-sigma uncertainty)

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
                i = int(k[1:])-1
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

total = len(good_plates)
#poop = np.zeros(5)
#indi = np.zeros(5)
#poopy = np.zeros(5)
#vdisp = np.zeros(total)
#vrot = np.zeros(total)
#disp_rot = np.zeros(5)


filename = 'outflow_galaxies'
with PdfPages(filename) as pdf:
    
    for i in range(0,445):
        plate, ifu = good_plates[i].split('-')
        name = mpl7_dir + 'HYB10-GAU-MILESHC/' + str(plate) + '/' + str(ifu) + '/manga-' + str(plate) + '-' + str(ifu) + '-MAPS-HYB10-GAU-MILESHC.fits.gz'
        if os.path.isfile(name):
            match = np.where((drpdata.field('plate') == int(plate)) & (drpdata.field('ifudsgn') == str(ifu)))
            hdu_dap = fits.open(name)
            print('('+str(i)+')', plate, ifu)
            #indi[i] = plate

            # hdu_dap.info()
            emlc = channel_dictionary(hdu_dap, 'EMLINE_SFLUX')
            dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Ha-6564'], :, :]
            dap_ha_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[emlc['Ha-6564'], :, :]
            dap_hb_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Hb-4862'], :, :]
            dap_ha_vel = hdu_dap['EMLINE_GVEL'].data[emlc['Ha-6564'], :, :]
            dap_ha_sig = hdu_dap['EMLINE_GSIGMA'].data[emlc['Ha-6564'], :, :]
            s_n_ha = dap_ha_sflux * np.sqrt(dap_ha_sivar)
        

            # read in the position angle and axis from from the NSA
            theta = drpdata.field('nsa_sersic_phi')[match][0]
            ba = drpdata.field('nsa_sersic_ba')[match][0]
            inc = np.arccos(ba)


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
                dap_serr = np.sqrt(1. / dap_sivar)
                dap_ha_serr = np.sqrt(1./dap_ha_sivar)

                # getting rid of bad spaxels
                # s_n = dap_sflux / dap_serr
                # bad = (s_n < 2.5)
                # s_n[bad] = 0
                # dap_vel[bad] = -2000
                # for b in range(0, len(dap_vel[i])):
                #     if dap_vel[b].all() < -1999:
                #         dap_vel.remove(dap_vel[b].all())
                #         dap_vel_ivar.remove(dap_vel[b].all())

                # creating vel_flip map
                vel_flip = np.zeros([size, size])
                ivar_flip = np.zeros([size, size])

                for m in range(0, size):
                    for n in range(0, size):
                        dist = np.sqrt((xproj_kpc - xproj_kpc[m, n]) ** 2 + (yproj_kpc + yproj_kpc[m, n]) ** 2)
                        test = (dist == np.min(dist))
                        vel_flip[m, n] = dap_vel[test]
                        ivar_flip[m, n] = dap_vel_ivar[test]

                delta_vel = dap_vel-vel_flip
                denom = np.sqrt((1/dap_ivar)+(1/ivar_flip))
                frac = delta_vel/denom
                flat = np.std(frac)

                # defining halfight radius + 1 arcsec
                re_arcsec = (drpdata.field('nsa_elpetro_th50_r')[match]) + np.array([1])
                re_kpc = re_arcsec * u.arcsec / cosmo.arcsec_per_kpc_proper(z)

                #Creating assymetry parameter
                re_large = np.ravel(radkpc) > re_kpc
                re_large_positive = np.ravel(yproj_kpc_map)[re_large] > 0
                re_large_negative = np.ravel(yproj_kpc_map)[re_large] < 0
                poop_positive = (np.std(np.ravel(frac)[re_large][re_large_positive]))
                poop_negative = (np.std(np.ravel(frac)[re_large][re_large_negative]))
                ha_good_large_positive = np.ravel(s_n_ha)[re_large][re_large_positive] > 5
                ha_good_large_negative = np.ravel(s_n_ha)[re_large][re_large_negative] > 5
                poop_positive_good = (np.std(np.ravel(frac)[re_large][re_large_positive][ha_good_large_positive]))
                poop_negative_good = (np.std(np.ravel(frac)[re_large][re_large_negative][ha_good_large_negative]))
                poop_tmp_good = (poop_positive_good + poop_negative_good) / 2.
                poop_tmp = (poop_positive + poop_negative) / 2.
                print('assymetry', poop_tmp)
                print('the assymetry parameter with good spaxels is', poop_tmp_good)
                #poop[i] = poop_tmp
                #poopy[i] = poop_tmp_good

                #Velocity dispersion for Vdisp/Vrot 
                ha_good_large = np.ravel(s_n_ha)[re_large] > 5
                med_vel_disp = np.median(np.ravel(dap_sig)[re_large][ha_good_large])
                print('the velocity dispersion is', med_vel_disp)
                #vdisp[i] = med_vel_disp

                #Velocity roation for Vdisp/Vrot
                major = yproj_kpc_map < 1.0 # within 1 kpc of major axis
                vel_major = (dap_vel)[major] # might need np.ravel() ?
                vrot = np.percentile(vel_major,90) # 90th percentile make sense?
                vrot_max = np.round(np.mean(vrot, dtype=None),5)
                print('the rotational velocity is', vrot_max)
                #vrot[i] = vrot_max

                #Vrot/Vdisp
                med_vel_disp_vrot_max = med_vel_disp/vrot_max
                print('the velocity dipsersion to rotation ratio is', med_vel_disp_vrot_max)
                #disp_rot[i] = med_vel_disp_vrot_max


                if (med_vel_disp_vrot_max > 0.3) and (poop_tmp_good > 1.8):
                    state = 'this galaxy has outflows!!'
                    print('this galaxy has outflows!')
                else:
                    state = 'no outflows detected'
                    print('no outflows detected')

                #plotting galaxies with outflows   
                if state == 'this galaxy has outflows!!':

                    #Identifying bad spaxels
                    s_n_ha = dap_ha_sflux/dap_ha_serr
                    bad = (s_n_ha < 3.)
                    s_n_ha[bad]=0

                    s_n = dap_sflux/dap_serr
                    bad = (s_n < 3.)
                    s_n[bad]=0

                    #flipping velocity and fining the asymmetry
                    vel_flip = np.zeros([size, size])
                    ha_hb = dap_ha_sflux / dap_hb_sflux
                    ha_hb_flip =np.zeros([size,size])
                    ha_vel_flip = np.zeros([size, size])

                    for m in range(0, size):
                        for n in range(0, size):
                            dist = np.sqrt((xproj_kpc - xproj_kpc[m, n]) ** 2 + (yproj_kpc + yproj_kpc[m, n]) ** 2)
                            test = (dist == np.min(dist))
                            vel_flip[m, n] = dap_vel[test]
                            ha_hb_flip[m, n] = ha_hb[test]
                            ha_vel_flip[m, n]= dap_ha_vel[test]

                    flip_difference = dap_vel-vel_flip
                    flip_ha_difference = dap_ha_vel-ha_vel_flip

                    flip_ha_difference[bad]= -2000
                    for b in range(0,size):
                        if flip_ha_difference[b].all() < -1999:
                            flip_ha_difference.remove(flip_ha_difference[b].all())

                    flip_difference[bad]= -2000
                    for b in range(0,size):
                        if flip_difference[b].all() < -1999:
                            flip_difference.remove(flip_difference[b ].all())

                    #Identifying outflow velocity
                    dap_vel_out = flip_difference/(2*np.cos(inc))
                    dap_vel_out_avg = np.median(dap_vel_out)
                    print(dap_vel_out_avg)

                    dap_vel_ha_out = flip_difference/(2*np.cos(inc))
                    dap_vel_ha_out_avg = np.median(dap_vel_out)
                    print(dap_vel_ha_out_avg)

                    # page 1, plot 1: galaxy coordinates in kpc
                    fig = plt.figure()
                    ax = fig.add_subplot(3,3,1)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(zproj_kpc_map, origin='lower', interpolation='nearest', cmap=cm.coolwarm)
                    plt.colorbar()
                    plt.title('vertical distance [kpc]', fontsize=10)
                    plt.suptitle('('+str(i)+')'+' '+plate+'-'+ifu + ', pg 1')

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
                    #ax = fig.add_subplot(3,3,6)
                    #ax.set_xlim(0, size)
                    #ax.set_ylim(0, size)
                    #ax.set_xticklabels(())
                    #ax.set_yticklabels(())
                    #plt.imshow(dap_sig, origin='lower',
                           #interpolation='nearest',
                           #cmap=cm.YlOrRd, vmin=0, vmax=250)
                    #plt.colorbar()
                    #plt.title('GSIGMA', fontsize=10)

                
                    # page 1, plot 7: emission-line velocity
                    ax = fig.add_subplot(3,3,7)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_vel, origin='lower', interpolation='nearest', cmap=cm.coolwarm, vmin=-250, vmax=250)
                    plt.colorbar()
                    plt.title('GVEL', fontsize=10)

                    # page 1, plot 7: Outflow velocity
                    ax = fig.add_subplot(3, 3, 6)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_vel_out, origin='lower', interpolation='nearest', cmap=cm.coolwarm, vmin=-250, vmax=250)
                    plt.colorbar()
                    plt.title('Outflow Velocity', fontsize=8)

                    # page 1, plot 8: Vflip
                    ax = fig.add_subplot(3,3,8)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(vel_flip, origin='lower', interpolation='nearest', cmap=cm.coolwarm, vmin=-250, vmax=250)
                    plt.colorbar()
                    plt.title('VFLIP', fontsize=10)

                    #page 1, plot 9: GVEL-Vflip
                    ax = fig.add_subplot(3,3,9)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(flip_difference, origin='lower', interpolation='nearest', cmap=cm.coolwarm, vmin=-80, vmax=80)
                    plt.colorbar()
                    plt.title('GVEL-VFLIP', fontsize=10)

    
                    pdf.savefig()
                    plt.close()

os.system("open %s &" % filename)







