from marvin.tools.maps import Maps
import os
import marvin.utils.plot.map as mapplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob
import numpy as np

#maps = Maps(plateifu='7977-12704', mode='remote')
#maps = Maps(plateifu='7977-12704', mode='local')
#maps = Maps(plateifu='7977-12704', mode='local', filename='/usr/netapp/physics/linux-lab/data/sas/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/7977/12704/manga-7977-12704-MAPS-SPX-GAU-MILESHC.fits.gz')
#maps = Maps(plateifu='7977-12704', mode='local', data_origin='file')
#print(maps)
# get an emission line map

sas_dir = os.environ['SAS_BASE_DIR']
plate = '7977'

maps_file = glob.glob(sas_dir+'/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/'+plate+'/*/manga-*MAPS-SPX-GAU*')


filename='ad_marvin_subplot.pdf'
with PdfPages(filename) as pdf:

    for i in range(0,1): #len(maps_file)):
        maps = Maps(filename=maps_file[i])
        print(maps)
      
        fig = plt.figure()

        ax = fig.add_subplot(2,2,1)
        haf = maps['emline_sflux_ha_6564']
        mapplot.plot(dapmap=haf, fig=fig, ax=ax)
        
        cb_kws={}
        cb_kws['vmin']=0.1
        cb_kws['vmax']=10.
        cb_kws['ticks']=[0.3,1,3]
        cb_kws['size']=10.
        ax = fig.add_subplot(2,2,2)
        o2f = maps['emline_sflux_oiid_3728']
        mapplot.plot(dapmap=o2f, fig=fig, ax=ax, cb_kws=cb_kws, log_cb=1, cbrange=(0.1,100))

        ax = fig.add_subplot(2,2,3)
        trash={}
        trash['size']=10
        trash['label']='Your MUM'
        junk={}
        junk['labelsize']=10
        cb_kws['label_kws']=trash
        cb_kws['tick_params_kws']=junk
        blah = getattr(o2f,'value')
        logo2f=np.log10(getattr(o2f,'value')*10)
        #mapplot.plot(value=logo2f, cbrange=(0,3), fig=fig, ax=ax, cb_kws=cb_kws)
        ticks=np.array([0.1,1,10,100])
        cb_kws['ticks']=ticks
        mapplot.plot(dapmap=o2f, fig=fig, ax=ax, cb_kws=cb_kws, log_cb=1, cbrange=(0.1,100))
        #plt.colorbar(vmin=0, vmax=100)

        #ax = fig.add_subplot(2,2,3)
        #o2v = maps['emline_gvel_oiid_3728']
        #mapplot.plot(dapmap=o2v, fig=fig, ax=ax) 
        #
        #ax = fig.add_subplot(2,2,4)
        #o2s = maps['emline_gsigma_oiid_3728']
        #mapplot.plot(dapmap=o2s, fig=fig, ax=ax, log_cb=1) 

        plt.suptitle(maps_file[i])
        fig.tight_layout()
        pdf.savefig()
        plt.close()

os.system("evince %s &" % filename)
