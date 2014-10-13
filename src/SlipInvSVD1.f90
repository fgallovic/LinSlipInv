!   Singular Value Decomposition applied to the linear slip inversion problem
!   AUTHOR:  Frantisek Gallovic
    
! WARNING: USE CULA ONLY IN COMBINATION WITH MKL !

    PROGRAM SlipInvSVD1
    USE SISVDmodule
#ifdef CULA
    USE cula_status
    USE cula_lapack
#endif
    IMPLICIT NONE
    REAL*8,ALLOCATABLE:: GTd(:),work(:),GTdvec(:,:,:),eigvec(:,:,:)
    INTEGER lwork,info
    INTEGER i,j,k
#ifdef CULA
    REAL*4,ALLOCATABLE:: culaG(:,:),culaW(:),culaVT(:,:)
#endif

    CALL Init()
    allocate(G(Nsvd,Msvd),D(Nsvd),W(Msvd))
    if(NRseis>0)allocate(normdat(NSTAcomp))
    if(NRgps>0)allocate(sigmaGPS(3,NRgps))
#ifdef MKL
    allocate(VT(minMNsvd,Msvd))
#else
    allocate(V(Msvd,Msvd))
#endif
    CALL CreateGandD()

    write(*,*)'  (saving matrix G and vector D to G.dat ...)'
    open(111,form='unformatted',FILE='G.dat')
    do i=1,Msvd
      write(111)G(1:Nsvd,i)
    enddo
    write(111)D
    close(111)

    write(*,*)'Calculating vector GTd ...'
    allocate(GTd(Msvd),GTdvec(Ssvd,NL,NW))
#ifdef MKL
    call dgemv('T',Nsvd,Msvd,1.d0,G ,Nsvd,D,1,0.d0,GTd,1)
#else
    GTd=matmul(D,G)
#endif
    GTdvec=RESHAPE(GTd,(/Ssvd, NL, NW/))
    write(*,*)'  (saving)'
    open(198,FILE='GTd.dat')
    write(198,'(1E13.5)')GTd
    close(198)
    open(198,FILE='GTd1d.dat')
    do i=1,Ssvd
      write(198,'(1000E13.5)')(sum(GTdvec(i,j,:)),j=1,NL)
    enddo
    close(198)
    deallocate(GTd,GTdvec)

    write(*,*)'Decompozition of matrix',Nsvd,'x',Msvd,'...'
#ifndef MKL
    write(*,*)'  (using Numerical Recepies)'
    CALL svdcmp(G,Nsvd,Msvd,Nsvd,Msvd,W,V) ! Warning: G became U
    goto 233
#endif
#ifdef CULA
    write(*,*)'  (using GPU...)'
    allocate(culaG(Nsvd,Msvd),culaW(Msvd),culaVT(minMNsvd,Msvd))
    culaG=G
    info = cula_initialize()
    if(info.ne.0)stop
    write(*,*) '  (...CULA initiated)'
    if(minSVDchoice==1)then    ! Warning: G became U
      info = cula_sgesvd('O','S',Nsvd,Msvd,culaG,Nsvd,culaW,culaG,Nsvd,culaVT,minMNsvd)
    else
      info = cula_sgesvd('O','A',Nsvd,Msvd,culaG,Nsvd,culaW,culaG,Nsvd,culaVT,minMNsvd)
    endif
    if(info.ne.0)then
      CALL cula_check_status(info)
      stop
    endif
    write(*,*) '  (CULA finished)'
    call cula_shutdown()

    do i=1,Msvd   !Assuming that the singular vectors are sorted descendently
      if(culaW(i)<1.e-10)exit      ! CULA WORKS IN REAL*4 AND THUS HAS PROBLEMS WITH VERY SMALL SINGULAR VALUES
    enddo
    culaW(i:Msvd)=0.
    
    G=culaG
    W=culaW
    VT=culaVT
    deallocate(culaG,culaW,culaVT)
    goto 233
#endif
#ifdef MKL
    write(*,*)'  (using MKL)'
    allocate(work(1))
    lwork=-1
	if(minSVDchoice==1)then
      CALL dgesvd('O','S',Nsvd,Msvd,G,Nsvd,W,G,Nsvd,VT,minMNsvd,work,lwork,info)
    else
      CALL dgesvd('O','A',Nsvd,Msvd,G,Nsvd,W,G,Nsvd,VT,minMNsvd,work,lwork,info)
    endif
    lwork=int(work(1))
    deallocate(work)
    allocate(work(lwork))
    if(minSVDchoice==1)then    ! Warning: G became U
      CALL dgesvd('O','S',Nsvd,Msvd,G,Nsvd,W,G,Nsvd,VT,minMNsvd,work,lwork,info)
    else
      CALL dgesvd('O','A',Nsvd,Msvd,G,Nsvd,W,G,Nsvd,VT,minMNsvd,work,lwork,info)
    endif
    if(info.ne.0)stop
    deallocate(work)
#endif
   
233 write(*,*)'Done. Saving ...'

    write(*,*)'  (singular values)'
    open(198,FILE='singularvalues.dat');write(198,'(1E12.5)')W(:);close(198)

    write(*,*)'  (SVD.dat)'
    open(111,form='unformatted',FILE='SVD.dat')
    write(111)D
    if(NRseis>0)write(111)normdat
    if(NRgps>0)then
      write(111)normdatGPS
      write(111)sigmaGPS
    endif
    write(111)W
#ifdef MKL
    do i=1,minMNsvd
      write(111)VT(i,1:Msvd)   ! MKL SVD returns VT instead of V
    enddo
#else
    do i=1,minMNsvd
      write(111)V(1:Msvd,i)
    enddo
#endif
    i=0;write(111)i
    write(111)G
    i=0;write(111)i
    close(111)

    END

