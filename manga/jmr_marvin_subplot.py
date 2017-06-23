# Jose Ruiz 20170621, goal is to add more plots
from marvin.tools.maps import Maps
import os
import marvin.utils.plot.map as mapplot
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


filename='jmr_marvin_subplot.pdf'
with PdfPages(filename) as pdf:

    for i in range(0, len(maps_file)):
        maps = Maps(filename=maps_file[i])
        print(maps)
        # page one plots multiple things, which we will plot multiple examples for each in the next pages
        fig = plt.figure()
        # h-alpha emission-gas
        ax = fig.add_subplot(2,2,1)
        haf = maps['emline_sflux_ha_6564']
        mapplot.plot(dapmap=haf, fig=fig, ax=ax)
        #[O II] emission-gas
        ax = fig.add_subplot(2,2,2)
        o2f = maps['emline_sflux_oiid_3728']
        mapplot.plot(dapmap=o2f, fig=fig, ax=ax)        
        # plot of [O II] velocity
        ax = fig.add_subplot(2,2,3)
        o2v = maps['emline_gvel_oiid_3728']
        mapplot.plot(dapmap=o2v, fig=fig, ax=ax) 
        # velocity dispersion [O II]
        ax = fig.add_subplot(2,2,4)
        o2s = maps['emline_gsigma_oiid_3728']
        mapplot.plot(dapmap=o2s, fig=fig, ax=ax) 

        plt.suptitle(maps_file[i])
        # the number corresponds to their wavelength in A*
        pdf.savefig()
        plt.close()
        # page 2
        #All of this figures represent gas emission
        fig=plt.figure()
        ax = fig.add_subplot(2,2,1)
        hbf = maps['emline_sflux_hb_4862']
        mapplot.plot(dapmap=hbf, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,2)
        o1f = maps['emline_sflux_oi_6302']
        mapplot.plot(dapmap=o1f, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,3)
        n2f= maps['emline_sflux_nii_6549']
        mapplot.plot(dapmap=n2f, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,4)
        s2f = maps['emline_sflux_sii_6718']
        mapplot.plot(dapmap=s2f, fig=fig, ax=ax)
        pdf.savefig()
        plt.close()
        #page 3
        # All the following figures represent velocities corresponding to each gas emission above
        fig=plt.figure()
        ax = fig.add_subplot(2,2,1)
        hbv = maps['emline_gvel_hb_4862']
        mapplot.plot(dapmap=hbv, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,2)
        o1v = maps['emline_gvel_oi_6302']
        mapplot.plot(dapmap=o1v, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,3)
        n2v = maps['emline_gvel_nii_6549']
        mapplot.plot(dapmap=n2v, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,4)
        s2v = maps['emline_gvel_sii_6718']
        mapplot.plot(dapmap=s2v, fig=fig, ax=ax)
        pdf.savefig()
        plt.close()
        #4 All the following figures represent the velocity dispersion
        fig=plt.figure()
        ax = fig.add_subplot(2,2,1)
        hbd = maps['emline_gsigma_hb_4862']
        mapplot.plot(dapmap=hbd, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,2)
        o1d = maps['emline_gsigma_oi_6302']
        mapplot.plot(dapmap=o1d, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,3)
        n2d = maps['emline_gsigma_nii_6549']
        mapplot.plot(dapmap=n2d, fig=fig, ax=ax)
        ax = fig.add_subplot(2,2,4)
        s2d = maps['emline_gsigma_sii_6718']
        mapplot.plot(dapmap=s2d, fig=fig, ax=ax)
        
        
        pdf.savefig()
        plt.close()

        
os.system("open %s &" % filename)
