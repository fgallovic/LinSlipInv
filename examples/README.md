#Example files for testing the inversion codes
----------------------------------------------

####List of available examples:
 - `SIV1a`: Benchmark from the Source Inversion Validation (SIV, see the enclosed pdf file, or directly http://equake-rc.info/sivdb/).
 
####Examples to become available soon:
 - `LAquila-realdata`: Inversion of real data recorded during the Mw6.3 2009 L'Aquila earthquake. You are welcome to try all the features for analysis of the slip resolution by means of synthetic tests.

####How to run the examples:
 1. Create a working directory a copy there all the files from the present directory.
 2. First step is the calculation of the Green's functions: Copy all files from the source directory src-dwn to the working directory and run first script  `firststep.sh` and then `calculate.sh`. You should eventually have all the Green's in file `NEZsor.dat`.
 3. Now you can start with the data inversion: Compile the source files in the directory `src` and create links to the executables in the working directory. Run the inversion (e.g., `SlipInvNNLS`).
 4. To plot the results, copy everything from the directory `src-graphics` to the working directory. Then you can run the specific plotting scripts and compare the plots with those in the directory `results.NNLS` in the specific example directory.
 5. You can start with experiments regarding different inversion approaches and synthetic/resolution tests.

Good luck!
