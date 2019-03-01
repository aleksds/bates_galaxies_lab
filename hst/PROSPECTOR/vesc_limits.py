import numpy as np
from astropy.constants import G
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)


mass_lo = np.array([1.1e9,7.9e9]) * const.M_sun
re_hi_arc = np.array([0.033,0.029]) * u.arcsec
z = np.array([0.603,0.459])

for i in range(0,len(z)):
    
    re_hi_kpc = re_hi_arc[i] / cosmo.arcsec_per_kpc_proper(z[i]) 
    vesc_lo = np.sqrt(G*mass_lo[i]/re_hi_kpc).to('km/s')
    print(vesc_lo)
