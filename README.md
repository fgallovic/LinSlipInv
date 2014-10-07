#LinSlipInv
===========

Linear multi time-window earthquake slip inversion with k^-2 smoothing

This code is a suite of codes for linear slip inversions and resolution analysis.

####Capabilities of the codes:
 - inversion of provided data for a given fault geometry
 - resolution analysis by means of synthetic tests with prescribed target rupture model or slip-rate pulse model
 - inversion in-depth analysis by means of spectral analysis of the forward matrix **G**.

####Implemented data types for inversions:
 - seismic waveforms including processed HR-GPS
 - static GPS vectors

####Implemented regularizations of the inversion:
 - truncated SVD
 - spatial *k*^-2 prior covariance function

####Included codes for evaluation of Green's functions:
 - Axitra (full-wavefield in 1D layered media)
 - Okada (static displacements in homogeneous halfspace)
