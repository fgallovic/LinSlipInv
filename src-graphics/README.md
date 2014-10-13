#Graphical output
-----------------

Fortran 90 codes that prepare Gnuplot scripts and its input files for plotting rupture models and seismograms.

####List of codes:
 - `mtilde2anime.f90`: Plots the inverted rupture evolution in terms of snapshots.
 - `slipratesonfault.f90`: Plots slip rate functions along the fault.
 - `SlipInv-seisplot.f90`: Plots seismogram comparison.

Notes:
 - Each of the codes is supplemented by its own script that compiles the code, executes it and calles `gnuplot` for plotting.
 - Only input files used for the inversion codes are need (no additional ones are required).
