#Inversion codes
----------------

Each of the inversion codes consists of several source files. Some of the source files are shared.

####Linear multi time-window earthquake slip inversion with *k*^-2 smoothing and positivity constraint (NNLS).
 - The code can be compiled by the following scripts `compile.SlipInvNNLS.sh`.

####Eigenanalysis and Truncated SVD solution of the linear multi time-window earthquake slip inversion
 - `SlipInvSVD1`: Applies SVD to the linear multi time-window earthquake slip inversion, providing singular values and vectors. Compile using script `compile.SlipInvSVD1.sh`.
 - `SlipInvSVD2`: Uses the singular values and vectors calculated by `SlipInvSVD1` to solve the slip inversion. Compile using script `compile.SlipInvSVD2.sh`

Notes:
 - Input files are shared by all the codes.
 - The use of MKL is higly recommended. SVD is several orders faster than that from the Numerical recipes.

Credits:
 - Intel MKL library for utlimate performance in linear algebra (SVD, matrix operations).
 - Nonnegative Least Square (NNLS) code originally developed by Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory (and published in the book "SOLVING LEAST SQUARES PROBLEMS", Prentice-Hall, 1974, translated to Fortran 90 by Alan Miller (February 1997).
 - Modified NNLS subroutine by Luo and Duraiswami (2011) taking advantage of OpenMP and the Intel MKL library.
 - Subroutine XAPIIR (IIR filter design and implementation) by Dave Harris (1990).
 - Subroutines Numerical recipes (FFT, SVD).
 - CULA subroutine for SVD on GPU.
 
 