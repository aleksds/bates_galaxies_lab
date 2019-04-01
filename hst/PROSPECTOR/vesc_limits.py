import numpy as np
from astropy.constants import G
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# should read this from a table
mass_lo = np.array([1.1e9,7.9e9,2.6e9,8.2e9,2.3e9,2.9e9,1.0e10,3.2e9,7.4e9,4.7e9,9.5e9,1.4e10]) * const.M_sun
re_hi_arc = np.array([0.033,0.029,0.018,0.023,0.041,0.058,0.021,0.023,0.088,0.177,0.044,0.045]) * u.arcsec
z = np.array([0.603,0.459,0.712,0.514,0.467,0.451,0.658,0.608,0.402,0.449,0.728,0.752])
galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']

for i in range(0,len(z)):
    
    re_hi_kpc = re_hi_arc[i] / cosmo.arcsec_per_kpc_proper(z[i]) 
    vesc_lo = np.sqrt(G*mass_lo[i]/re_hi_kpc).to('km/s')
    print(galaxies[i], re_hi_kpc, vesc_lo)
