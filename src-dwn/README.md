#Axitra
-------

Code based on discrete wavenumber method providing full-wavefield Green's functions in 1D layered media 
(Kennett and Kerry, 1979; Bouchon, 1981; Coutant, 1989).

The code was originally written by O. Coutant.
Additional modifications have been made by J. Zahradník, J. Burjánek and F. Gallovič.

###How to use the code

AXITRA consists of four codes. The best approach is to use batch files `firststep.sh` and `calculate.sh`.
The first file compiles all the codes and runs code `prepare.f90` for preparation of the AXITRA calculations,
in particular it prepares list of elementary sources covering the rupture in regular grid.
Then, a parallel loop using `xargs` is to be started by `calculate.sh`.
The number of processors can be set in the batch file.
For each elementary source the codes `gr_nez.for` and `cnv_nez.for` are run automatically.
Intermediate results including Green’s functions (GFs) for the individual elementary
sources are stored in the dat directory. Finally, the Green’s functions are resorted
by `resort.f90` into the `NEZsor.dat` file. Note that the order of the GFs is such
that the outer loop is over stations. Thus if more crustal models are to be considered,
the `NEZsor.dat` files for the individual subsets of stations and then simply appended one after each other.

###References:
- Kennett, B. L. N., and N. J. Kerry (1979). Seismic waves in a stratified half
space, Geophys. J. Roy. Astron. Soc. 57, 557–583.
- Bouchon, M. (1981). A simple method to calculate Green’s functions for
elastic layered media, Bull. Seismol. Soc. Am. 71, 959–971.
- Coutant, O. (1989). Program of Numerical Simulation AXITRA, Research
Report, Lab. de Geophys. Interne et Tectonophys., Grenoble, France
