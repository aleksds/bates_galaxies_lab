#This code graphs the r_e from 875 divided by the r_e from 475 vs the r_e from 875 to illustrate whether or not there's a systematic difference between r_e at the two wavelengths. I used the values found from the coarse images in both filters except using the fine psf's.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os


#re in order:J0826,J0901,J0905,J0944,J1107,J1219,J1341,J1506,J1558,J1613,J2116,J2140
four = [0.25,0.6,0.24,0.51,0.21,0.57,0.59,0.22,0.64,2.08,0.34,0.25]
eight = [0.4,0.51,0.23,0.24,0.61,0.83,0.28,0.31,3.07,3.29,0.62,0.52]
eightoverfour=np.zeros(12)
for i in range(12):
    eightoverfour[i] = eight[i]/four[i]

filename = 'cl_wavelength_re_comp.pdf'
with PdfPages(filename) as pdf:
    fig = plt.figure()
    plt.title('re_F814W / re_F475W vs re_F814W' ,fontsize=20)
    plt.ylabel('re_F814W/re_F475W', fontsize=17)
    plt.xlabel('re_F814W (pix)', fontsize=17)

    markers = ['x','D','o','>','.',',','v','1','2','8','*','+']
    for c in range(12):
        plt.scatter(eight[c],eightoverfour[c], marker = markers[i],s=30)
    a=plt.scatter(eight[0],eightoverfour[0], marker = markers[0],s=30)
    b=plt.scatter(eight[1],eightoverfour[1], marker = markers[1],s=30)
    c=plt.scatter(eight[2],eightoverfour[2], marker = markers[2],s=30)
    d=plt.scatter(eight[3],eightoverfour[3], marker = markers[3],s=30)
    e=plt.scatter(eight[4],eightoverfour[4], marker = markers[4],s=30)
    f=plt.scatter(eight[5],eightoverfour[5], marker = markers[5],s=30)
    g=plt.scatter(eight[6],eightoverfour[6], marker = markers[6],s=30)
    h=plt.scatter(eight[7],eightoverfour[7], marker = markers[7],s=30)
    i=plt.scatter(eight[8],eightoverfour[8], marker = markers[8],s=30)
    j=plt.scatter(eight[9],eightoverfour[9], marker = markers[9],s=30)
    k=plt.scatter(eight[10],eightoverfour[10], marker = markers[10],s=30)
    l=plt.scatter(eight[11],eightoverfour[11], marker = markers[11],s=30)
    plt.legend((a,b,c,d,e,f,g,h,i,j,k,l),
           ('J0826','J0901','J0905','J0944','J1107','J1219','J1341','J1506','J1558','J1613','J2116','J2140'),
           scatterpoints=1,
           loc='upper left',
           ncol=2,
           fontsize=8)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
os.system('open %s &' % filename)
