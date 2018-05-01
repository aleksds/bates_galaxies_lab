#This is a program to plot the magnitiude of different models for the J0905 galaxy

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#sersic values:
#These values are from independant trials for each filter
re_F475W_indep_s = 0.365
re_F814W_indep_s = 0.212

#These values are from galfitm where x and y values are allowed to be different between filters
re_F475W_simul_s = 1.388*10**(-17)
re_F814W_simul_s = 0.09853  

#These values are from galfitm where x and y positions must remain constant
re_F475W_simul2_s = 1.110*10**(-16)
re_F814W_simul2_s = 0.108

re_diff_s = [(re_F475W_indep_s/re_F814W_indep_s), (re_F475W_simul_s/re_F814W_simul_s), (re_F475W_simul2_s/re_F814W_simul2_s)]
F814W_s = [re_F814W_indep_s, re_F814W_simul_s, re_F814W_simul2_s]

with PdfPages('ec_re_plot.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0905 Difference in Effective Radius Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('re_F475W/reF814W')

    labels = ['sersic independent', 'sersic simultaneous', 'sersic semi-simultaneous', 'psf independent', 'psf simultaneous', 'psf semi-simultaneous']
    colors = ['lightcoral', 'firebrick', 'darkorange', 'palegreen', 'lightseagreen', 'lightskyblue']
    for i in range(0,3):
        plt.scatter(F814W_s[i],re_diff_s[i], label = labels[i], color = colors[i])
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.25))


    pdf.savefig(dpi=1000)
    plt.close()

os.system('open %s &' % 'ec_re_plot.pdf')
