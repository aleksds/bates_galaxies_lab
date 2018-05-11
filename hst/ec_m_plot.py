#This is a program to plot the magnitiude of different models for the J0905 galaxy

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#sersic values:
#These values are from independant trials for each filter
m_F475W_indep_s = 19.497
m_F814W_indep_s = 19.089

#(semi-simultaneous) These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_s = 19.537  
m_F814W_simul_s = 19.101

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_s = 19.396
m_F814W_simul2_s = 18.966

m_diff_s = [(m_F475W_indep_s - m_F814W_indep_s), (m_F475W_simul_s - m_F814W_simul_s), (m_F475W_simul2_s - m_F814W_simul2_s)]
F814W_s = [m_F814W_indep_s, m_F814W_simul_s, m_F814W_simul2_s]

#psf values:
#These values are from independant trials for each filter
m_F475W_indep_p = 19.626
m_F814W_indep_p = 19.224

#(semi-simultaneous) These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_p = 19.597  
m_F814W_simul_p = 19.203

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_p = 19.663
m_F814W_simul2_p = 19.231

m_diff_p = [(m_F475W_indep_p - m_F814W_indep_p), (m_F475W_simul_p - m_F814W_simul_p), (m_F475W_simul2_p - m_F814W_simul2_p)]
F814W_p = [m_F814W_indep_p, m_F814W_simul_p, m_F814W_simul2_p]

with PdfPages('ec_m_plot.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0905 Difference in Magnitude Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('m_F475W - m_F814W')
    labels = ['sersic independent', 'sersic semi-simultaneous', 'sersic simultaneous', 'psf independent', 'psf semi-simultaneous', 'psf simultaneous']
    colors = ['lightcoral', 'firebrick', 'darkorange', 'palegreen', 'lightseagreen', 'lightskyblue']
    markers = ['o', 'v','s','o', 'v','s']
    for i in range(0,6):
        if i < 3:
            plt.scatter(F814W_s[i],m_diff_s[i], label = labels[i], color = colors[i], marker = markers[i])
        if i >= 3:
            plt.scatter(F814W_p[i-3],m_diff_p[i-3], label = labels[i], color = colors[i], marker = markers[i])
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.25))


    pdf.savefig(dpi=1000)
    plt.close()

os.system('open %s &' % 'ec_m_plot.pdf')

#on plot, values go from left to right as follows: indep, simul, simul2

