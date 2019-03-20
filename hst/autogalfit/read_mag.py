# Aleks Diamond-Stanic, March 20, 2019
# goal: read information from galfit output files to determine color and magnitude information
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

dir = sys.argv[1]
gal = sys.argv[2]
files = glob.glob(dir+'/'+gal+'*.band')
print(files)
chi = np.zeros(len(files))
umag = np.zeros(len(files))
vmag = np.zeros(len(files))
jmag = np.zeros(len(files))
sizepix = np.zeros(len(files))

for i in range(0, len(files)):
    print(files[i])
    with open(files[i]) as f:
        content = f.readlines()
        print(content[48])
        print(content[47])
        chi[i] = np.float(content[3][14:19])
        vmag[i] = np.float(content[47][4:10])
        umag[i] = np.float(content[47][11:17])
        jmag[i] = np.float(content[47][18:24])
        sizepix[i] = np.float(content[48][4:9])
        print(chi[i], umag[i], vmag[i], jmag[i], sizepix[i])

print('F475W: ', np.min(umag), np.mean(umag), np.median(umag), np.max(umag), np.std(umag), (np.max(umag)-np.min(umag))/2.)

print('F814W: ', np.min(vmag), np.mean(vmag), np.median(vmag), np.max(vmag), np.std(vmag), (np.max(vmag)-np.min(vmag))/2.)

print('F160W: ', np.min(jmag), np.mean(jmag), np.median(jmag), np.max(jmag), np.std(jmag), (np.max(jmag)-np.min(jmag))/2.)

        

name = 'mag_color_'+gal+'.pdf'

with PdfPages(name) as pdf:   
    fig = plt.figure()

    plt.scatter(sizepix, umag, marker='o', color='blue')
    plt.scatter(sizepix, vmag, marker='+', color='green')
    plt.scatter(sizepix, jmag, marker='s', color='red')
    plt.ylim([20.5, 17.5])
    plt.xlim([0, np.median(sizepix)*3])
    plt.ylabel('Magnitude')
    plt.xlabel('Half-light radius in pixels')
    plt.title(dir+'/'+gal)
    
    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.scatter(vmag-jmag, umag-vmag, marker='o', color='black')
    plt.xlabel('[F814W] - [F160W]')
    plt.ylabel('[F475W] - [F814W]')
    plt.xlim([-0.5, 1.0])
    plt.ylim([-0.5, 1.0])

    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.scatter(sizepix, umag-vmag, marker='o', color='black')
    plt.xlabel('Half-light radius in pixels')
    plt.ylabel('[F475W] - [F814W]')
    plt.xlim([0, np.median(sizepix)*3])
    plt.ylim([-0.5, 1.0])    

    pdf.savefig()
    plt.close()

    fig = plt.figure()

    plt.scatter(sizepix, vmag-jmag, marker='o', color='black')
    plt.xlabel('Half-light radius in pixels')
    plt.ylabel('[F814W] - [F160W]')
    plt.xlim([0, np.median(sizepix)*3])
    plt.ylim([-0.5, 1.0])    

    pdf.savefig()
    plt.close()
    
os.system('open %s &' % name)
    
