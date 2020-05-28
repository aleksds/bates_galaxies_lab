
# coding: utf-8

# In[8]:


import os
import numpy as np
import copy

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from marvin.tools import Maps
from marvin import config
import marvin.utils.plot.map as mapplot
from marvin.tools.image import Image

# config.access = 'collab'
# config.mode = 'remote'
# config.login()
# config.setRelease('MPL-8')
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)


# In[9]:


def snr_mask(f_map, snr_cut):
    flux_mask = (np.abs(f_map.value * np.sqrt(f_map.ivar)) < snr_cut)
    return flux_mask


# In[12]:


def flip_map(map_obj, original_map, flux_map, minimum_snr): 
    flipped_flux = copy.deepcopy(flux_map)
    flipped_map  = copy.deepcopy(original_map)
    asymm_map    = copy.deepcopy(original_map)
    
    # Retrieving NSA values
    ba, theta, z,  re_arcesc = (map_obj.nsa.get(key) for key in ["sersic_ba", "sersic_phi", "z", "elpetro_th50_r"])
    inc = np.arccos(ba)
    re_arcsec = re_arcesc + np.array([1])
    
    #create arrays that correspond to x and y coordinates for each spaxel
    size = len(original_map)
    major_axis = np.zeros([size, size], dtype=bool) 
            
    xpos = np.zeros(original_map.shape); 
    ypos = np.zeros(original_map.shape)
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
    radkpc = radpix / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

    # compute values along the x and y axis of the disk and the z axis above the disk
    xproj_kpc = xproj / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
    yproj_kpc = yproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
    zproj_kpc = zproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

    radkpc_map = np.zeros([size, size])
    xproj_kpc_map = np.zeros([size, size])
    yproj_kpc_map = np.zeros([size, size])
    zproj_kpc_map = np.zeros([size, size])

    for j in range(0, size):
        for k in range(0, size):
            radkpc_map[j, k] = radkpc[j, k] / u.kpc
            yproj_kpc_map[j, k] = yproj_kpc[j, k] / u.kpc
            zproj_kpc_map[j, k] = zproj_kpc[j, k] / u.kpc

    for m in range(0, size):
        for n in range(0, size):
            dist = np.sqrt((xproj_kpc - xproj_kpc[m, n]) ** 2 + (yproj_kpc + yproj_kpc[m, n]) ** 2)
            test = (dist == np.min(dist))

            if (np.ravel(np.where(test == True)) == [m,n]).all():
                major_axis[m, n] = True

            if size == len(original_map):
                flipped_map.value[m, n] = original_map.value[test]
                flipped_map.ivar[m, n] = original_map.ivar[test]
                flipped_map.mask[m, n] = original_map.mask[test]

                flipped_flux.value[m, n] = flux_map.value[test]
                flipped_flux.ivar[m, n] = flux_map.ivar[test]
            else:
                pass
    
    mask_flipped_flux = snr_mask(flipped_flux, minimum_snr)
    
    asymm_map.value[:,:] = original_map - flipped_map
    asymm_map.value[(original_map.value ==0) & (flipped_map.value ==0)] = 0 
    asymm_map.ivar = 1. / ((1. / original_map.ivar) + (1. / flipped_map.ivar))
    asymm_map.mask = mask_flipped_flux | snr_mask(flux_map, minimum_snr)

    return flipped_map, mask_flipped_flux, asymm_map, major_axis


# In[13]:


def disk(map_obj, original_map, disk_re_size=1.):
    # Retrieving NSA values
    ba, theta, z,  re_arcesc = (map_obj.nsa.get(key) for key in ["sersic_ba", "sersic_phi", "z", "elpetro_th50_r"])
    inc = np.arccos(ba)
    re_arcsec = re_arcesc + np.array([1])
     
    #create arrays that correspond to x and y coordinates for each spaxel
    size = len(original_map)
            
    xpos = np.zeros(original_map.shape); 
    ypos = np.zeros(original_map.shape)
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
    radkpc = radpix / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

    # compute values along the x and y axis of the disk and the z axis above the disk
    xproj_kpc = xproj / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
    yproj_kpc = yproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
    zproj_kpc = zproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

    radkpc_map = np.zeros([size, size])
    xproj_kpc_map = np.zeros([size, size])
    yproj_kpc_map = np.zeros([size, size])
    zproj_kpc_map = np.zeros([size, size])

    for j in range(0, size):
        for k in range(0, size):
            radkpc_map[j, k] = radkpc[j, k] / u.kpc
    
    re_kpc = (disk_re_size * re_arcsec) * u.arcsec / cosmo.arcsec_per_kpc_proper(z)        
    cen = abs(radkpc_map) < (re_kpc / u.kpc)
    return cen

