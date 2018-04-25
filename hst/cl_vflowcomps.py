#this code makes a plot of the observed escape velocities vs the toy model (using heckman 2011 equation) escape velocities
#it also makes a plot of the observed escape velocities vs the escape velocities we calculated
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

vflow = [1228,1206,2470,1778,1828,1830,875,0,829,2416,1456,606]
toy = [[7440,5890],[3319,3597],[7991,8134],[4343,6347],[6915,4070],[5155,4270],[7005,10149],[11586,9783],[2717,1240],[3523,2801],[4963,3682],[7332,5088]]
esc = [[1758,1392],[1337,1481],[1148,1169],[1489,2176],[2193,1291],[2540,2104],[993,1439],[1377,1163],[1955,892],[1854,1474],[2026,1503],[1103,766]]

markers = ['*','x']
colors = ['blue','green']

x = np.linspace(0,3000,3000)
y = np.linspace(0,3000,3000)

with PdfPages('cl_vflowvstoy.pdf') as pdf:
    plt.figure()
    plt.plot(x,y,linestyle=':',color='red',linewidth=1, label='x = y')
    plt.xlim(0,2500)
    plt.ylim(0,12000)
    for w in range (0,12):
        for i in range(0,2):
            if w == 0:
                if i == 0:
                    plt.scatter(vflow[w], toy[w][0], s=30, marker = markers[0], c=colors[0], label='F475W')
                if i == 1:
                    plt.scatter(vflow[w], toy[w][1], s=30, marker = markers[1], c=colors[1], label='F814W')
            else:
                if i == 0:
                    plt.scatter(vflow[w], toy[w][i], s=30, marker = markers[0], c=colors[0])
                if i == 1:
                    plt.scatter(vflow[w], toy[w][i], s=30, marker = markers[1], c=colors[1])
                
    plt.legend(loc='upper right', prop={'size': 12})
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Toy Model Speed (km/s)')
    plt.title('Observed Outflow Speed vs. Toy Model Speed')
    pdf.savefig()
    plt.close()

markers = ['v','o']
# now for esc vel

#list of escape velocities in this order: best, low, high (also first is 475 then 814)
vels = [[[1758,1413,2265],[1392,1119,1794]],[[1337,1111,1662],[1481,1204,1802]],[[1148,871,1446],[1169,887,1471]],[[1489,1183,1832],[2176,1728,2677]],[[2193,1783,2729],[1291,1049,1606]],[[2540,2089,3125],[2104,1730,2589]],[[2540,2089,3125],[2104,1730,2589]],[[993,798,1236],[1439,1156,1791]],[[1377,1057,1714],[1163,892,1447]],[[1955,1571,2461],[892,717,1123]],[[1854,1525,2229],[1474,1212,1772]],[[2026,1628,2229],[1503,1208,1982]],[[1103,857,1506],[766,594,1045]]]

with PdfPages('cl_vflowvsescvel.pdf') as pdf:
    plt.figure()
    plt.xlim(500,3000)
    plt.ylim(-0.1,6.3)
    #plt.plot(x,y,linestyle=':',color='red',linewidth=1)
    for w in range (0,12):
        for i in range(0,2):
            
            if w == 0:
                if i == 0:
                    plt.scatter(vflow[w], vflow[w]/vels[w][0][0], s=30, marker = markers[0], c=colors[0], label = 'F475W')
                    error = [np.array([vflow[w]/(vels[w][0][2])]),np.array([vflow[w]/(vels[w][0][1])])]
                    plt.errorbar(vflow[w], vflow[w]/vels[w][0][0], yerr = error,linewidth=0.5,color=colors[0])
                else:
                    plt.scatter(vflow[w], vflow[w]/vels[w][1][0], s=30, marker = markers[1], c=colors[1], label = 'F814W')
                    error = [np.array([vflow[w]/(vels[w][1][2])]),np.array([vflow[w]/(vels[w][1][1])])]
                    plt.errorbar(vflow[w], vflow[w]/vels[w][1][0], yerr = error,linewidth=0.5,color=colors[1])
                
            else:
                plt.scatter(vflow[w], vflow[w]/vels[w][i][0], s=30, marker = markers[i], c=colors[i])
                error = [np.array([vflow[w]/(vels[w][i][2])]),np.array([vflow[w]/(vels[w][i][1])])]
                plt.errorbar(vflow[w], vflow[w]/vels[w][i][0], yerr = error,linewidth=0.5,color=colors[i])
           

                    
                    
    plt.legend(loc='upper left', prop={'size': 12})
    plt.xlabel('Observed Outflow Speed (km/s)')
    plt.ylabel('Observed Outflow Speed / Estimated Escape Velocity')
    plt.title('Observed Outflow Speed vs. Estimated Escape Velocity')
    pdf.savefig()
    plt.close()

    
