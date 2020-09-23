import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

galaxies =              ['J0826','J0901','J0905','J0944','J1107','J1219','J1341','J1506','J1558','J1613','J2116','J2140']

vavg = np.array([         1131.2, 1210.1, 2395.3, 1172.4, 1359.8, 1543.7,  766.9, 1378.9, 807.1, 1437.5, 1056.2,  342.8])
vmax = np.array([         1456.2, 1575.4, 2884.8, 1878.3, 2015.6, 1962.6, 1995.2, 1768.8,1151.2, 2374.1, 1915.1,  950.4])

v50_davis = np.array(   [1247.52,1308.73,2514.93,1286.62,1416.76,1684.03, 760.78,1288.30,868.45, 807.53,1124.12, 512.42])
v95_davis = np.array(   [1657.28,1777.88,3019.30,1988.41,2050.23,2031.88,1997.66,2447.93,1255.89,2503.18,2238.64,1042.26])


filename = 'vavg_vmax_tremonti_davis.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    vel_array = np.arange(101)*40

    plt.scatter(vmax, v95_davis-vmax, marker='o', color='#ff7f0e', label=r'$v_{max}$')
    plt.scatter(vavg, v50_davis-vavg, marker='+', color='#ff7f0e', label=r'$v_{avg}$')

    plt.ylim(-750,750)
    plt.xlim(0,3200)
    plt.ylabel(r'$|v_{Davis}|-|v_{Tremonti}|$ [km s$^{-1}$]', fontsize=13)
    plt.xlabel(r'Tremonti Outflow Velocity [km s$^{-1}$]', fontsize=13)
    plt.plot(vel_array, np.zeros(len(vel_array)), label=r'$v_{Tremonti} = v_{Davis}$', linestyle='dashed', color='#2ca02c')

    for i in range(0,len(vavg)):
        plt.text(vmax[i], v95_davis[i]-vmax[i], galaxies[i], fontsize=7)
        plt.text(vavg[i], v50_davis[i]-vavg[i],galaxies[i], fontsize=7)
    
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()

    print(np.median(v95_davis-vmax))
    print(np.median(v50_davis-vavg))
    
    plt.legend(fontsize=11, loc='lower right')

    pdf.savefig()
    plt.close()

    
os.system('open %s &' % filename)


