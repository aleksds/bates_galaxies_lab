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
files = glob.glob(dir+'/'+gal+'_'+band+'_output.galfit.*.band')
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

name = 'chi_values.pdf'

with PdfPages(name) as pdf:   
    fig = plt.figure()

    plt.scatter(mag, sizepix, marker='o', color='orange')
    plt.yscale('log')
    plt.ylim([0.1, 3.])
    plt.title(dir+'/'+gal+'_'+band)

    for i in range(0, len(files)):
        plt.text(mag[i], sizepix[i], str(chi[i]), fontsize=3)

    pdf.savefig()
    plt.close()

os.system('open %s &' % name)
    
