This archive contains the Python program voronoi_2d_binning.py
implementing the two-dimensional adaptive spatial binning method
of Cappellari & Copin (2003, MNRAS, 342, 345).

Usage example is provided by the procedure
voronoi_2d_binning_example.py.

Perform the following simple steps to bin you own 2D data with minimal
Python interaction:

1) Write your data vectors [X,Y,Signal,Noise] in the text file
voronoi_2d_binning_example.txt, following the example provided;

2) Change the line "targetSN = 50.0" in the procedure
voronoi_2d_binning_example.py, to specify the desired target S/N of
your final bins;

3) Run the program "voronoi_2d_binning_example" and wait for the final plot 
to appear. The output is saved in the text file voronoi_2d_binning_output.txt. 
The last column BIN_NUM is all is needed to actually bin the data.

4) Read the documentation at the beginning of the file
voronoi_2d_binning.py to fully understand the meaning of the various
optional output parameters.

----------------------------------
When some pixels have very low S/N
----------------------------------

2D-Binning should not be used blindly when some pixels contain
significant noise but virtually no signal. This situation may happen
e.g. when extracting the gas kinematics from observed galaxy spectra.
One way of using voronoi_2d_binning consists of first selecting the
pixels with S/N above a minimum threshold and then 2D-Binning each set
of connected pixels *separately*. Alternatively one may optimally weight 
the pixels before binning. See Sec.2.1 of Cappellari & Copin (2003) for
details.

---------------------
2D-Binning X-ray data
---------------------

For X-ray data, or other data coming from photon-counting devices the
noise is generally accurately Poissonian. In the Poissonian case the S/N
in a bin can never decrease by adding a pixel (see Sec.2.1 of Cappellari
& Copin 2003), and it is preferable to bin the data *without* first
removing the observed pixels with zero or very low signal.

--------------------------
2D-Binning very big images
--------------------------

Computation time in voronoi_2d_binning scales nearly as npixels^1.5, so
it may become a problem for large images (e.g. at the time of writing
npixels > 1000x1000). Let's assume that we really need to bin the image
as a whole, and that the S/N in a significant number of pixels is well
above our target S/N. As for many other computational problems, a way to
radically decrease the computation time consists of proceeding in a
hierarchical manner. Suppose for example we have a 4000x4000 pixels
image, we can do the following: (i) we rebin the image regularly (e.g.
in groups of 8x8 pixels) to a manageable size of 500x500 pixels; (ii) we
apply the standard Voronoi 2D-binning procedure to the 500x500 image;
(iii) we transform all unbinned pixels (which already had enough S/N) of
the 500x500 Voronoi 2D-binned image back into their original individual
full-resolution pixels; (iv) we now apply Voronoi 2D-binning only the
connected regions of full-resolution pixels; (v) we merge the set of
lower resolution bins with the higher resolution ones.


Michele Cappellari
cappellari_at_astro.ox.ac.uk
Vicenza, 13 February 2003

Last software update 31 March 2016
