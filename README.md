# bates_galaxies_lab

This repository includes code developed by members of the Bates
Galaxies Lab (sometimes called "the bagel," which in principle could
stand for Bates Astrophysics Galaxy Evolution Lab) at Bates College in
Lewiston, Maine.  Current lab members (as of the Fall 2016 semester),
include Aleks Diamond-Stanic (Assistant Professor of Physics and
Astronomy), Sophia Gottlieb (senior physics major, class of 2017), and
Joshua Rines (senior physics major, class of 2017).

# ad_centroid.py

This code compares centroid values based on three routines from
photutils: centroid_com, centroid_1dg, centroid_2dg.  The current
version of the code plots these centroid values on top of a postage
stamp image for the galaxy J0905.

# ad_compilation.py

This is a version of the jr_compilation*py code that calculates
fluxes, luminosities, colors, mass-to-light ratios, and stellar
masses.  This is an outdated piece of code that was used to check the
mass-to-light ratios and total stellar mass estimate for the galaxy
J0905.

# ad_phot_sed.py

This code performs photometry on three images
using circular apertures and then plots the spectral energy
distribution for each aperture.  The current version uses 13 different
apertures for the galaxy J0826 and includes an SED for each annulus.  

# bgl_aper_phot.py

This code makes "postage stamp" visualizations of a galaxy in three
different filters, performs aperture photometry on each image, and
plots the flux (energy/area/time) in circular apertures centered on
the galaxy.  The current version is set up to work for the galaxy
J0826.

# bgl_gau_bkg.py

This code analyzes the distribution of pixel values for an entire
image and quantifies the mean background flux and its dispersion by
fitting a Gaussian function to a binned histogram of pixel values.
The current version analyzes the F814W image for the galaxy J0826.

# bgl_image_stamp.py

This code reads in images of the same galaxy in three filters and
makes "postage stamp" visualizations with the goal of providing a
qualitative sense for how the morphology of the galaxy and its
luminosity change as a function of wavelength.  This code is expanded
upon in bgl_aper_phot.py and bgl_pix_vis.py.  The current version is
set up to work for the galaxy J0826.

# bgl_pix_vis.py

The code makes "postage stamp" visualizations of a galaxy in three
different filters and then produces color-magntidue diagrams that show
the relationship between flux in a given filter and the color measured
with respect to an adjacent filter for all pixels in the postage
stamps.  The current version is set up to work for the galaxy J0905.

# jr_*.py

Aweseome code being developed by Josh Rines

# jr_aper_phot.py

# jr_color_color.py

# jr_compilation.py

# jr_compilation_J0826.py

# jr_compilation_J0905.py

# jr_compilation_J1107.py

# jr_flux_plots.py

# jr_flux_plots_dp.py

# jr_flux_vs_wavelength.py

# jr_gau_bkg.py

# jr_image_gau.py

# jr_image_stamp.py

# jr_overlays_J0826.py

# jr_overlays_J0905.py

# jr_overlays_J1107.py

# jr_phot_sed.py

# jr_phot_sed_noloop.py

# jr_phot_sed_noloop_figure3.py

# jr_phot_sed_noloop_figure4.py

# sg_*.py

Awesome code being developed by Sophia Gottlieb

# sg_comp_overlay.py

# sg_flux_pix.py

# sg_gau_bkg.py

# sg_image_gau.py

# sg_image_stamp.py

# sg_phot_sed.py






