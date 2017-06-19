from marvin.tools.maps import Maps
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

maps = Maps(plateifu='7977-12704', mode='remote')
#maps = Maps(plateifu='7977-12704', mode='local')
print(maps)
# get an emission line map

filename='ad_marvin_test.pdf'
with PdfPages(filename) as pdf:

    fig = plt.figure()
    haflux = maps.getMap('emline_gflux', channel='ha_6564')
    values = haflux.value
    ivar = haflux.ivar
    mask = haflux.mask
    haflux.plot()
    
    pdf.savefig()
    plt.close()

os.system("evince %s &" % filename)
