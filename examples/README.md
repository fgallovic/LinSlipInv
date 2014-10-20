#Example files for testing the inversion codes
----------------------------------------------

####List of available examples:
 - `SIV1a`: Benchmark from the Source Inversion Validation (SIV, see the enclosed pdf file, or directly http://equake-rc.info/sivdb/).
 - `LAquila-realdata`: Inversion of real data recorded during the Mw6.3 2009 L'Aquila earthquake. You are welcome to try all the features for analysis of the slip resolution by means of synthetic tests.

####How to run the examples:
 1. Create a working directory a copy there all the files from the present directory.
 2. First step is the calculation of the Green's functions: Copy all files from the source directory src-dwn to the working directory and run first script  `firststep.sh` and then `calculate.sh`. You should eventually have all the Green's in file `NEZsor.dat`. In case of the real data, you can use the precomputed
    3D GFs considering 3D tomographic velocity model of Di Stefano et al. (2011) with real topography located in subdirectory `3D_GFs_topo_3D_GPS'.
 3. In case of the real data application, run conversion of real data files in subdirectory `data` in the example directory using script `processseis.sh`. Copy the resulting seismogram files `rvseis[nez].dat` or create their symbolic links into the working directory.
 3. Now you can start with the data inversion: Compile the source files in the directory `src` and create symbolic links to the executables in the working directory. Run the inversion (e.g., `SlipInvNNLS`).
 4. To plot the results, copy everything from the directory `src-graphics` to the working directory. Then you can run the specific plotting scripts and compare the plots with those in subdirectory `results.NNLS` in the example directory.
 5. You can start with experiments regarding different inversion approaches and synthetic/resolution tests. The latter can be attained by modyfing input file
    `SlipInvSVD.in`.

Good luck!

####Credits:
 - Observed waveforms were provided by the Italian Strong Motion Database (ITACA, http://itaca.mi.ingv.it/). 
   The AQU station belongs to the MedNet network (http://mednet.rm.ingv.it/data.php),
   while the remaining ones belong to the Italian strong-motion network (RAN)
   (http://www.protezionecivile.gov.it/jcms/it/ran.wp).
 - The displacement waveforms obtained from HR-GPS stations were kindly provided by A. Avallone
   (Avallone et al., 2011).
 - 3D GFs were calculated for the 3D tomographic velocity model of Di Stefano et al. (2011) with
   real topography (Gallovič et al., 2014).

####References
 - Avallone, A., M. Marzario, A. Cirella, A. Piatanesi, A. Rovelli, C. Di Alessandro,
   E. D’Anastasio, N. D’Agostino, R. Giuliani, and M. Mattone (2011). Very high rate (10 Hz)
   GPS seismology for moderate-magnitude earthquakes: The case of the M w 6.3 L’Aquila (central Italy)
   event, J. Geophys. Res., 116(B2), B02305, doi:10.1029/2010JB007834.
 - Di Stefano, R., C. Chiarabba, L. Chiaraluce, M. Cocco, P. De Gori, D. Piccinini, and L. Valoroso (2011).
   Fault zone properties affecting the rupture evolution of the 2009 (Mw 6.1)
   L’Aquila earthquake (central Italy): insights from seismic tomography, Geophys. Res. Lett. 38, L10310.
 - Gallovič, F., Imperatori, W., Mai, P. M. (2014). Effect of three-dimensional velocity heterogeneities
   and topography on slip inversions: case study of the Mw6.3 2009 Lâ€™Aquila earthquake,
   submitted to J. Geophys. Res.

