#LinSlipInv
-----------

Linear multi time-window earthquake slip inversion with *k*<sup>-2</sup> smoothing

Suite of codes for linear slip inversions and resolution analysis.

####Capabilities of the codes:
 - Inversion of provided data for a given (possibly segmented) fault geometry
 - Resolution analysis by means of synthetic tests with prescribed target rupture model or slip-rate pulse model
 - Inversion in-depth analysis by means of spectral analysis of the forward matrix **G**
 - Can take advantage of Intel MKL library and/or CULA (GPU) for faster performance
 - Plotting the results

####Possible data types for inversions:
 - Seismic waveforms including processed HR-GPS
 - Static GPS vectors

####Implemented regularizations of the inversion:
 - Truncated SVD
 - Spatial *k*<sup>-2</sup> prior covariance function
 - Positivity constraint on the slip rate functions

####Included codes for evaluation of Green's functions:
 - Axitra (full-wavefield in 1D layered media)
 - Okada (static displacements in homogeneous halfspace)

------------

###Content of directories:
 - `src` - Inversion codes
 - `src-stations` - Converts stations locations from lat,long to X,Y (X points towards north, Y towards east)
 - `src-dwn` - Axitra code for Green's function calculations
 - `src-graphics` - Codes for generating graphics (requires Gnuplot)
 - `examples` - Several examples for testing the code
 - `papers` - Papers related to the inversion codes, explaining basics of the SVD and NNLS approaches, resolution analysis, etc.
 