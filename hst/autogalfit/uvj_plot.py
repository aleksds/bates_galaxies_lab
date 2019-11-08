import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

extdata = ascii.read('bgl_uvj_ext.dat')
nucdata = ascii.read('bgl_uvj_nuc.dat')
totdata = ascii.read('bgl_uvj_tot.dat')

filename = 'uvj_plot.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    #test_u_flux = 10.**(totdata['umag'])
    
    ext_uv_color = extdata['umag']-extdata['vmag']
    ext_vj_color = extdata['vmag']-extdata['jmag']

    nuc_uv_color = nucdata['umag']-nucdata['vmag']
    nuc_vj_color = nucdata['vmag']-nucdata['jmag']

    tot_uv_color = totdata['umag']-totdata['vmag']
    tot_vj_color = totdata['vmag']-totdata['jmag']
    
    ext_uv_err = np.sqrt(extdata['uerr']**2+extdata['verr']**2)
    ext_vj_err = np.sqrt(extdata['verr']**2+extdata['jerr']**2)

    nuc_uv_err = np.sqrt(nucdata['uerr']**2+nucdata['verr']**2)
    nuc_vj_err = np.sqrt(nucdata['verr']**2+nucdata['jerr']**2)

    tot_uv_err = np.sqrt(totdata['uerr']**2+totdata['verr']**2+0.03**2)
    tot_vj_err = np.sqrt(totdata['verr']**2+totdata['jerr']**2+0.03**2)
    
    #plt.scatter(ext_vj_color, ext_uv_color, color='red')
    plt.errorbar(ext_vj_color, ext_uv_color, xerr=ext_vj_err, yerr=ext_uv_err, fmt='o', marker=',', elinewidth=1, color='red', label='extended')

    #plt.scatter(nuc_vj_color, nuc_uv_color, color='blue', marker='*')
    plt.errorbar(nuc_vj_color, nuc_uv_color, xerr=nuc_vj_err, yerr=nuc_uv_err, fmt='o', marker='*', elinewidth=1, color='blue', label='central')

    #plt.scatter(tot_vj_color, tot_uv_color, color='green', marker='o')
    plt.errorbar(tot_vj_color, tot_uv_color, xerr=tot_vj_err, yerr=tot_uv_err, fmt='o', marker='.', elinewidth=1, color='green', label='total')
    
    
    plt.xlabel('V-J')
    plt.ylabel('U-V')
    
    plt.xlim([-0.5,2.0])
    plt.ylim([0,2.5])

    plt.legend()
    
    # 0.807 - 1.998
    slant_vj = np.arange(1600-807+1)/1000+0.807
    slant_uv = 0.88 * slant_vj + 0.59
    plt.plot(slant_vj, slant_uv, color='black')

    horiz_vj = np.arange(500.+807+1)/1000.-0.5
    horiz_uv = np.zeros(len(horiz_vj))+1.3
    plt.plot(horiz_vj, horiz_uv, color='black')

    vert_uv = np.arange(2500-1998+1)/1000+1.998
    vert_vj = np.zeros(len(vert_uv))+1.6
    plt.plot(vert_vj, vert_uv, color='black')

    
    #for i in range(0, len(extdata)):
    #    plt.text(ext_vj_color[i], ext_uv_color[i], extdata['Galaxy'][i])
    #    plt.text(nuc_vj_color[i], nuc_uv_color[i], nucdata['Galaxy'][i])
    #    plt.text(tot_vj_color[i], tot_uv_color[i], totdata['Galaxy'][i])
    
    pdf.savefig()
    plt.close()

    

os.system('open %s &' % filename)
