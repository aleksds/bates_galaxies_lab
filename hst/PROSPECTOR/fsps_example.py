import os
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.2)

wave10, spec10 = sp.get_spectrum(tage=10)
wave3, spec3 = sp.get_spectrum(tage=3)
wave1, spec1 = sp.get_spectrum(tage=1)
wave03, spec03 = sp.get_spectrum(tage=0.3)
wave01, spec01 = sp.get_spectrum(tage=0.1)
wave003, spec003 = sp.get_spectrum(tage=0.03)
wave001, spec001 = sp.get_spectrum(tage=0.01)

filename = 'fsps_example.pdf'

with PdfPages(filename) as pdf:
    fig = plt.figure()
    plt.plot(wave10,spec10,color='red', label='10 Gyr')
    plt.plot(wave3,spec3,color='orange', label='3 Gyr')
    plt.plot(wave1,spec1,color='cyan', label='1 Gyr')
    plt.plot(wave03,spec03,color='green', label='300 Myr')
    plt.plot(wave01,spec01,color='blue', label='100 Myr')
    plt.plot(wave003,spec003,color='indigo', label='30 Myr')
    plt.plot(wave001,spec001,color='violet', label='10 Myr')
    plt.xlim(3000,10000)
    plt.ylim(1e-18,1e-12)
    plt.xscale("log", nonposx='clip')
    plt.yscale("log", nonposy='clip')
    plt.legend(loc='lower right', prop={'size': 10})
    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)
