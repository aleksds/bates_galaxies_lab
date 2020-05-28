import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter

r = (np.arange(10000)+1) / 100.

re = 1.

surf = np.exp(-7.669*((r/re)**0.25-1))

flux = np.zeros(len(r))
flux[0] = surf[0] * 2 * np.pi * r[0] * r[0]
for i in range(1,len(r)):
    flux[i] = flux[i-1] + surf[i] * 2. * np.pi * r[i] * (r[i] - r[i-1])
#flux = surf * (np.pi * r**2)

flux = flux / np.max(flux)
    
filename = 'devauc.pdf'

with PdfPages(filename) as pdf:
    
    fig = plt.figure()

    ax = plt.subplot(111)
    
    plt.scatter(r, surf)

    for axis in [ax.xaxis, ax.yaxis]:
      axis.set_major_formatter(ScalarFormatter())
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-9,1e3)
    
    pdf.savefig()
    plt.close()

    fig = plt.figure()

    ax = plt.subplot(111)
    
    plt.scatter(r, flux)

    for axis in [ax.xaxis, ax.yaxis]:
      axis.set_major_formatter(ScalarFormatter())
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-3,1e2)
    
    pdf.savefig()
    plt.close()

    
os.system('open %s &' % filename)

re_kpc = np.array([0.10118411, 0.08678413, 0.07553863, 0.06131937, 0.09174529, 0.14820475, 0.08846379, 0.07937191, 0.20856418, 0.74144918, 0.15673792, 0.10650737])

factor = 1. / re_kpc
for i in range(0,len(re_kpc)):
    element = int(round(factor[i]*100))
    print(flux[element])
    #print(flux[element]/0.5)

# fraction of light within r=1kpc
frac = np.array([0.960717691783886, 0.9706289023555622, 0.9778115612718672, 0.9858813411790087, 0.9673131546882592, 0.925344843140101, 0.9695044975464989, 0.9754351856431778, 0.8781548653578193, 0.5862226539058183, 0.9186408452041525, 0.9569267737212989])

corr = frac / 0.5
print('min: ', np.min(corr))
print('max: ', np.max(corr))
print('median: ', np.median(corr))
print('mean: ', np.mean(corr))
print('std: ', np.std(corr))



