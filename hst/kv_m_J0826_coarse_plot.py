#This is a program to plot the magnitiude of different models for the J0826 galaxy

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#sersic values:
#These values are from independant trials for each filter
m_F475W_indep_s = 19.486
m_F814W_indep_s = 18.992

#These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_s = 19.474 
m_F814W_simul_s = 19.037

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_s = 19.824
m_F814W_simul2_s = 18.824

m_diff_s = [(m_F475W_indep_s - m_F814W_indep_s), (m_F475W_simul_s - m_F814W_simul_s), (m_F475W_simul2_s - m_F814W_simul2_s)]
F814W_s = [m_F814W_indep_s, m_F814W_simul_s, m_F814W_simul2_s]

#psf values:
#These values are from independant trials for each filter
m_F475W_indep_p = 19.686
m_F814W_indep_p = 19.253

#These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_p = 20.236 
m_F814W_simul_p = 19.880

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_p = 20.283
m_F814W_simul2_p = 19.858

m_diff_p = [(m_F475W_indep_p - m_F814W_indep_p), (m_F475W_simul_p - m_F814W_simul_p), (m_F475W_simul2_p - m_F814W_simul2_p)]
F814W_p = [m_F814W_indep_p, m_F814W_simul_p, m_F814W_simul2_p]

with PdfPages('kv_m_J0826_coarse_plot.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0826 Difference in Magnitude Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('m_F475W - m_F814W')
##    plt.scatter(F814W_s, m_diff_s, label = 'sersic')
##    plt.scatter(F814W_p, m_diff_p, label = 'psf')
##    plt.legend(loc=3)
    labels = ['sersic independent (r_e, x, y, magnitude allowed to vary)', 'sersic simultaneous (x, y, magnitude allowed to vary)', 'sersic semi-simultaneous (magnitude allowed to vary)', 'psf independent', 'psf simultaneous', 'psf semi-simultaneous']
    colors = ['lightcoral', 'firebrick', 'darkorange', 'palegreen', 'lightseagreen', 'lightskyblue']
    marker = ['o', 's', 'D', 'o', 's', 'D']
    for i in range(0,6):
        if i < 3:
            plt.scatter(F814W_s[i],m_diff_s[i],s = 60, label = labels[i], color = colors[i], marker = marker[i])
        if i >= 3:
            plt.scatter(F814W_p[i-3],m_diff_p[i-3],s = 60, label = labels[i], color = colors[i], marker = marker[i])
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.25))


    pdf.savefig(dpi=1000)
    plt.close()

os.system('open %s &' % 'kv_m_J0826_coarse_plot.pdf')

#on plot, values go from left to right as follows: indep, simul, simul2
