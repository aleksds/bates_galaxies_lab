from marvin.tools.maps import Maps
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob

#maps = Maps(plateifu='7977-12704', mode='remote')
#maps = Maps(plateifu='7977-12704', mode='local')
#maps = Maps(plateifu='7977-12704', mode='local', filename='/usr/netapp/physics/linux-lab/data/sas/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/7977/12704/manga-7977-12704-MAPS-SPX-GAU-MILESHC.fits.gz')
#maps = Maps(plateifu='7977-12704', mode='local', data_origin='file')
#print(maps)
# get an emission line map

sas_dir = os.environ['SAS_BASE_DIR']
plate = '7977'

maps_file = glob.glob(sas_dir+'/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/'+plate+'/*/manga-*MAPS-SPX-GAU*')


filename='ad_marvin_test.pdf'
with PdfPages(filename) as pdf:

    for i in range(0, len(maps_file)):
        maps = Maps(filename=maps_file[i])
        print(maps)
      
        fig = plt.figure()
    
        haflux = maps.getMap('emline_gflux', channel='ha_6564')
        values = haflux.value
        ivar = haflux.ivar
        mask = haflux.mask
        haflux.plot()
        plt.suptitle(maps_file[i])
        
        pdf.savefig()
        plt.close()

os.system("open %s &" % filename)
