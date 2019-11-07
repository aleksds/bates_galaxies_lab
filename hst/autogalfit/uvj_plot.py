import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

data = ascii.read('bgl_uvj_ext.dat')

filename = 'uvj_plot.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    uv_color = data['umag']-data['vmag']
    vj_color = data['vmag']-data['jmag']

    uv_err = np.sqrt(data['uerr']**2+data['verr']**2)
    vj_err = np.sqrt(data['verr']**2+data['jerr']**2)
    
    plt.scatter(vj_color, uv_color)
    plt.errorbar(vj_color, uv_color, xerr=vj_err, yerr=uv_err, fmt='o', elinewidth=1)

    plt.xlabel('V-J')
    plt.ylabel('U-V')
    
    plt.xlim([0,2.5])
    plt.ylim([0,2.5])

    # 0.807 - 1.998
    slant_vj = np.arange(1600-807+1)/1000+0.807
    slant_uv = 0.88 * slant_vj + 0.59
    plt.plot(slant_vj, slant_uv, color='black')

    horiz_vj = np.arange(807+1)/1000.
    horiz_uv = np.zeros(len(horiz_vj))+1.3
    plt.plot(horiz_vj, horiz_uv, color='black')

    vert_uv = np.arange(2500-1998+1)/1000+1.998
    vert_vj = np.zeros(len(vert_uv))+1.6
    plt.plot(vert_vj, vert_uv, color='black')

    
    #for i in range(0, len(data)):
    #    plt.text(vj_color[i], uv_color[i], data['Galaxy'][i])
    
    pdf.savefig()
    plt.close()

    

os.system('open %s &' % filename)
