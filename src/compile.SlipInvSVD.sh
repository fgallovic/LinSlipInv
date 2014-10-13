#SlipInvSVD1
#-----------
#Compile with MKL and CUDA (GPU)
#ifort -fpp -mkl -DMKL -DCULA $CULA_INC_PATH/cula_status.f90 $CULA_INC_PATH/cula_lapack.f90 -oSlipInvSVD1CULA SlipInvSVD1.f90 CreateGandD.f90 filters.for dc3dmodif.f init.f90 nr.for $linkCULA

#Compile with MKL
#ifort -O -fpp -mkl -DMKL -oSlipInvSVD1 SlipInvSVD1.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for
ifort -O -fpp -mkl -DMKL -oSlipInvSVD2 SlipInvSVD2.f90 init.f90 OutputModel.f90

#Compile without MKL support
ifort -O -fpp -oSlipInvSVD1 SlipInvSVD1.f90 CreateGandD.f90 filters.for init.f90 dc3dmodif.f nr.for
#ifort -O -fpp -oSlipInvSVD2 SlipInvSVD2.f90 init.f90 OutputModel.f90


#SlipInvSVD2
#-----------
#Compile with MKL

#Compile without MKL support
