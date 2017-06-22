# Cbradna, goal is to get ratio graph
from marvin.tools.maps import Maps
import os
import marvin.utils.plot.map as mapplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob

# maps = Maps(plateifu='7977-12704', mode='remote')
# maps = Maps(plateifu='7977-12704', mode='local')
# maps = Maps(plateifu='7977-12704', mode='local', filename='/usr/netapp/physics/linux-lab/data/sas/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/7977/12704/manga-7977-12704-MAPS-SPX-GAU-MILESHC.fits.gz')
# maps = Maps(plateifu='7977-12704', mode='local', data_origin='file')
# print(maps)
# get an emission line map

sas_dir = os.environ['SAS_BASE_DIR']
plate = '7977'

maps_file = glob.glob(
    sas_dir + '/mangawork/manga/spectro/analysis/v2_0_1/SPX-GAU-MILESHC/' + plate + '/*/manga-*MAPS-SPX-GAU*')

filename = 'jmr_marvin_subplot.pdf'
with PdfPages(filename) as pdf:
    for i in range(0, len(maps_file)):
        maps = Maps(filename=maps_file[i])
        print(maps)

        fig = plt.figure()
        # h-alpha emission-gas
        ax = fig.add_subplot(2, 2, 1)
        haf = maps['emline_sflux_ha_6564']
        mapplot.plot(dapmap=haf, fig=fig, ax=ax)
        # [O II] emission-gas
        ax = fig.add_subplot(2, 2, 2)
        o2f = maps['emline_sflux_oiid_3728']
        mapplot.plot(dapmap=o2f, fig=fig, ax=ax)
        # plot of [O II] velocity
        ax = fig.add_subplot(2, 2, 3)
        o2v = maps['emline_gvel_oiid_3728']
        mapplot.plot(dapmap=o2v, fig=fig, ax=ax)
        # velocity dispersion [O II]
        ax = fig.add_subplot(2, 2, 4)
        o2s = maps['emline_gsigma_oiid_3728']
        mapplot.plot(dapmap=o2s, fig=fig, ax=ax)

        plt.suptitle(maps_file[i])
        # the number corresponds to their wavelength in A*
        pdf.savefig()
        plt.close()
        # page 2
        fig = plt.figure()
        ax = fig.add_subplot(2, 2, 1)
        h_o_ratio = Maps.getMapRatio(maps,'emline_sflux','ha_6564','oiid_3728')
        # h_o_ratio = maps['ha_6564/oiid_3728']
        mapplot.plot(dapmap=h_o_ratio, fig=fig, ax=ax)

        pdf.savefig()
        plt.close()
os.system("open %s &" % filename)
