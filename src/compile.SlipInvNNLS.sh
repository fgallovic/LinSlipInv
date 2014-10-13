#Compile with MKL NNLS
icc -O -c -openmp nnlsmkl.c
ifort -O -fpp -openmp -mkl -DMKL -DNNLSMKL -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnlsmkl.o OutputModel.f90

#Compile with classical NNLS
#ifort -O -fpp -openmp -mkl -DMKL -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnls.f90 OutputModel.f90

#Compile without MKL support
#ifort -O -fpp -openmp -oSlipInvNNLS SlipInvNNLS.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for nnls.f90 OutputModel.f90
