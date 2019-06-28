from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import numpy as np

# J082638.41+430529.3 0.037
# J090133.42+031412.4 0.057
# J090523.59+575912.4 0.031
# J094417.84+093019.3 0.026
# J110702.87+041702.7 0.076
# J121955.77+033615.9 0.021
# J134136.79-032125.2 0.037
# J150636.29+540220.7 0.016
# J155811.23+395720.8 0.014
# J161332.52+283414.7 0.041
# J211625.14-063444.8 0.147
# J214000.49+120914.5 0.143

# From Table 6 (using R_V=3.1) from Schlafly et al. (2011)
# WFC3 F475W: Ab / E(B-V) = 3.248
# WFC3 F814W: Ab / E(B-V) = 1.536
# WFC3 F160W: Ab / E(B-V) = 0.512

# From Table 2 (column 3) of Yuan et al. (2013)
# FUV: 4.37
# NUV: 7.06

# From Table 2 (column 5) of Yuan et al. (2013)
# FUV: 6.892
# NUV: 6.738

# from page 16 of Bianchi 2013 quotes the following
# FUV: 8.06
# NUV: 7.95

ra =  ['08h26m38.41s', '09h01m33.42s', '09h05m23.59s', '09h44m17.84s', '11h07m02.87s', '12h19m55.77s', '13h41m36.79s', '15h06m36.29s', '15h58m11.23s', '16h13m32.52s', '21h16m25.14s', '21h40m00.49s']
dec = ['43d05m29.3s',  '03d14m12.4s',  '57d59m12.4s',  '09d30m19.3s',  '04d17m02.7s',  '03d36m15.9s',  '-03d21m25.2s', '54d02m20.7s',  '39d57m20.8s',  '28d34m14.7s',  '-06d34m44.8s', '12d09m14.5s']
ebv = np.zeros(len(ra))
for i in range(0,len(ra)):
    coords = SkyCoord(ra[i], dec[i], frame='icrs')
    sfd = SFDQuery()
    ebv[i] = sfd(coords)
    print(ra[i], ': E(B-V) = {:.3f} mag'.format(ebv[i]), ', A_475 = {:.3f} mag'.format(ebv[i]*3.248), ', A_814 = {:.3f} mag'.format(ebv[i]*1.536), ', A_160 = {:.3f} mag'.format(ebv[i]*0.512))
    print('A_FUV = {:.3f} mag'.format(ebv[i]*4.37), ', A_NUV = {:.3f} mag'.format(ebv[i]*7.06) )
    print('A_FUV = {:.3f} mag'.format(ebv[i]*6.892), ', A_NUV = {:.3f} mag'.format(ebv[i]*6.738) )
    print('A_FUV = {:.3f} mag'.format(ebv[i]*8.06), ', A_NUV = {:.3f} mag'.format(ebv[i]*7.95) )

