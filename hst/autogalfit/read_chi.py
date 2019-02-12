# Aleks Diamond-Stanic, February 7, 2019
# goal: read information from galfit output files to construct chi-squared contours in terms of re and magnitude
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

dir = sys.argv[1]
gal = sys.argv[2]
band = sys.argv[3]
files = glob.glob(dir+'/'+gal+'_'+band+'*.band')
print(files)
chi = np.zeros(len(files))
mag = np.zeros(len(files))
sizepix = np.zeros(len(files))

for i in range(0, len(files)):
    print(files[i])
    with open(files[i]) as f:
        content = f.readlines()
        chi[i] = np.float(content[3][14:19])
        mag[i] = np.float(content[47][4:10])
        sizepix[i] = np.float(content[48][4:8])
        print(chi[i], mag[i], sizepix[i])

mag_1d = np.unique(mag)
sizepix_1d = np.unique(sizepix)
chi_2d = np.zeros([len(mag_1d), len(sizepix_1d)])
#count = 0
for i in range(0, len(mag_1d)):
    for j in range(0, len(sizepix_1d)):
        test = np.where((mag == mag_1d[i]) & (sizepix == sizepix_1d[j]))
        print(mag[test], sizepix[test], chi[test])
        chi_2d[j,i] = chi[test[0]]
        #print(mag[count], mag_1d[i])
        #print(sizepix[count], sizepix_1d[j]) 
        #count = count+1
        
        
        
name = 'chi_values_'+gal+'_'+band+'.pdf'

with PdfPages(name) as pdf:   
    fig = plt.figure()

    plt.scatter(mag, sizepix, marker='o', color='orange')
    plt.yscale('log')
    plt.ylim([0.1, 3.])
    plt.title(dir+'/'+gal+'_'+band)

    for i in range(0, len(files)):
        plt.text(mag[i], sizepix[i], str(chi[i]), fontsize=5)

    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.contourf(mag_1d, sizepix_1d, chi_2d, 20, cmap='RdGy')
    #plt.clim(0,20)
    plt.colorbar()
    
    pdf.savefig()
    plt.close()

    fig = plt.figure()

    #plt.contour(mag_1d, sizepix_1d, chi_2d)

    chi_min = np.min(chi_2d)
    #levels = np.array([chi_min, chi_min+1, chi_min+2.71, chi_min+4.00, chi_min+6.63, chi_min+9.00])##np.arange(10)/2 + chi_min
    levels = np.array([1.0, 4.00, 9.00])+chi_min
    cs = plt.contour(mag_1d, sizepix_1d, chi_2d, levels)
    plt.clabel(cs, inline=1, fontsize=10)
    plt.ylim([0, 2])
    plt.ylabel('Half-Light radius in pixels')
    plt.xlabel('Magnitude')
    labels = ['68%', '95%', '99.7%']
    #labels = ['68%', '90%','95%','99%','99.7%']
    for i in range(len(labels)):
        cs.collections[i].set_label(labels[i])

    plt.legend(loc='upper left')
    
    pdf.savefig()
    plt.close()
    
os.system('open %s &' % name)
    
