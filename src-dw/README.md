#Axitra
-------

Code based on discrete wavenumber method providing full-wavefield Green's functions in 1D layered media 
(Kennett and Kerry, 1979; Bouchon, 1981; Coutant, 1989).

The code was originally written by O. Coutant.
Additional modifications have been made by J. Zahradník, J. Burjánek and F. Gallovič.

###How to run the code
```
./firststep.sh
./calculate.sh
```
The first script compiles the required codes and runs `prepare.f90`.
The second script executes parallel (if possible) evaluation of Green's functions.
The result should be all Green's functions gathered in file `NEZsor.dat` required by the inversion codes.

###References:
- Kennett, B. L. N., and N. J. Kerry (1979). Seismic waves in a stratified half
space, Geophys. J. Roy. Astron. Soc. 57, 557–583.
- Bouchon, M. (1981). A simple method to calculate Green’s functions for
elastic layered media, Bull. Seismol. Soc. Am. 71, 959–971.
- Coutant, O. (1989). Program of Numerical Simulation AXITRA, Research
Report, Lab. de Geophys. Interne et Tectonophys., Grenoble, France
