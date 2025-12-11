#Compile with MKL NNLS
source load_intel
icx -O -c -qopenmp nnlsmkl.c
ifx -O -fpp -qopenmp -qmkl -DMKL -DNNLSMKL -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnlsmkl.o OutputModel.f90

#Compile with classical NNLS
#ifort -O -fpp -openmp -mkl -DMKL -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnls.f90 OutputModel.f90

#Compile without MKL support
#ifort -O -fpp -openmp -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnls.f90 OutputModel.f90
