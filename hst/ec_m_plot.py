#This is a program to plot the magnitiude of different models for the J0905 galaxy

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

m_F475W_indep = 20.492
m_F814W_indep = 20.469

#These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul = 19.337  
m_F814W_simul = 18.896

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2 = -0.221
m_F814W_simul2 = 19.125

m_diff = [(m_F475W_indep - m_F814W_indep), (m_F475W_simul - m_F814W_simul)]
F814W = [m_F814W_indep, m_F814W_simul]

with PdfPages('ec_m_plot.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0905 Difference in Magnitude Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.title('F475W')
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('m_F475W - m_F814W')
    plt.scatter(F814W, m_diff)

    pdf.savefig(dpi=1000)
    plt.close()

os.system('open %s &' % 'ec_m_plot.pdf')
