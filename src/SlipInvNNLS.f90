!   AUTHOR:  Frantisek Gallovic
    
    PROGRAM SlipInvNNLS
    USE SISVDmodule
#ifdef NNLSMKL
USE iso_c_binding
#endif
    IMPLICIT NONE
    REAL*8 rnorm
    REAL*8,ALLOCATABLE:: work(:)
    INTEGER,ALLOCATABLE:: iwork(:)
    INTEGER dumi,i

#ifdef NNLSMKL
! --- INTEROPERABLE VARIABLES
INTEGER,PARAMETER :: DP           = c_double ! ifdef USE_DOUBLE: c_double, otherwise: c_float
INTEGER(c_int) :: method          = 1        ! 0-1
INTEGER(c_int) :: Nnnls                      ! equations
INTEGER(c_int) :: Mnnls                      ! unknowns
INTEGER(c_int) :: NSYS            = 1        ! pocet soustav
INTEGER(c_int) :: isTransposed    = 0        ! matrix transpose
REAL(DP)       :: TOL_TERMINATION = 1d-6     ! tolerance
INTEGER(c_int) :: MKLThreads      = 6        ! MKL threads
INTEGER(c_int) :: OMPThreads      = 1        ! OpenMP threads
INTEGER(c_int) :: MAX_ITER_LS,MAX_ITER_NNLS
! --- INTERFACES TO C FUNCTIONS
INTERFACE

SUBROUTINE nnlsOMPSysMKL(A,b,x,isTransposed,maxNNLSIters,maxLSIters,nSys,m,n,MKLT,OMPT,TOL_TERMINATION) BIND(C,NAME='nnlsOMPSysMKL')
IMPORT c_int,DP
REAL(DP) :: A(*),b(*),x(*)
INTEGER(c_int),VALUE :: isTransposed,maxNNLSIters,maxLSIters,nSys,m,n,MKLT,OMPT
REAL(DP),VALUE :: TOL_TERMINATION
END SUBROUTINE

SUBROUTINE nnlsOMPSysMKLUpdates(A,b,x,isTransposed,maxNNLSIters,maxLSIters,nSys,m,n,MKLT,OMPT,TOL_TERMINATION) BIND(C,NAME='nnlsOMPSysMKLUpdates')
IMPORT c_int,DP
REAL(DP) :: A(*),b(*),x(*)
INTEGER(c_int),VALUE :: isTransposed,maxNNLSIters,maxLSIters,nSys,m,n,MKLT,OMPT
REAL(DP),VALUE :: TOL_TERMINATION
END SUBROUTINE

END INTERFACE
#endif


    CALL Init()
    allocate(G(Nsvd,Msvd),D(Nsvd),W(Msvd))
    if(NRseis>0)allocate(normdat(NSTAcomp))
    if(NRgps>0)allocate(sigmaGPS(3,NRgps))
    lambdanum=1
    allocate(Dout(Nsvd,lambdanum),M(Msvd,lambdanum))
    CALL CreateGandD()

    write(*,*)'Saving matrix G and vector D to G.dat ...'
    open(111,form='unformatted',FILE='G.dat')
    do i=1,Msvd
      write(111)G(1:Nsvd,i)
    enddo
    write(111)D
    close(111)

    write(*,*)'Calculating m_tilde using NNLS ...'
    write(*,*)'  (matrix ',Nsvd,'x',Msvd,')'

#ifdef NNLSMKL
    
    Nnnls=Nsvd
    Mnnls=Msvd
    MAX_ITER_LS=(Mnnls+Nnnls)*2
    MAX_ITER_NNLS=MAX_ITER_LS
    call nnlsOMPSysMKLUpdates(G,D,M(:,1),isTransposed,MAX_ITER_NNLS,MAX_ITER_LS,NSYS,Nnnls,Mnnls,MKLThreads,OMPThreads,TOL_TERMINATION)

#else

    allocate(work(Msvd),iwork(Msvd))
    CALL NNLS (G,Nsvd,Msvd,D,M(:,1),rnorm,work,iwork,dumi) 
    deallocate(work,iwork)

#endif

    open(111,form='unformatted',FILE='G.dat')
    do i=1,Msvd
      read(111)G(1:Nsvd,i)
    enddo
    read(111)D
    close(111)
    Dout=matmul(G,M)

    maxw=-1.d0;W=1.d0;lambdafrom=1 !ONLY FORMAL definitions

    CALL OutputModel()

    END
