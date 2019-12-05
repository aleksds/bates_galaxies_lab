import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import numpy as np


fraction = np.array([22/60, 17/59, 31/59, 26/59, 34/59, 34/59, 26/58, 34/59, 37/59, 42/60])
unc = np.array([np.sqrt(22)/60, np.sqrt(17)/59, np.sqrt(31)/59, np.sqrt(26)/59, np.sqrt(34)/59, np.sqrt(34)/59, np.sqrt(26)/58,
                    np.sqrt(34)/59, np.sqrt(37)/59, np.sqrt(42)/60])
sigma_sfr_p = np.log10(np.array([1.4e-4,7.6e-4,0.0013,0.0019,0.0025,0.0031,0.0041,0.0056,0.0087,0.041])) 

filename = 'ad_frac.pdf'

with PdfPages(filename) as pdf:

    fig = plt.figure()
    plt.errorbar(x=sigma_sfr_p, y=fraction, yerr=unc,ls='none')
   
    plt.scatter(sigma_sfr_p, fraction)
    plt.xlabel('Log(Surface Density of Star Formation)')
    plt.ylabel('Fraction of Edge-on Galaxies with Kinematic Asymmetry')
    plt.xlim(-4,-1)
    plt.ylim(0,1)
    z = np.polyfit(sigma_sfr_p, fraction, 1)
    p = np.poly1d(z)
    plt.plot(sigma_sfr_p,p(sigma_sfr_p),"r--")

    pdf.savefig()
    plt.close


os.system('open %s &' % filename)
