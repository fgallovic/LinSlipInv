!   Singular Value Decomposition applied to the linear slip inversion problem
!   AUTHOR:  Frantisek Gallovic
    
    PROGRAM SlipInvSVD2
    USE SISVDmodule
    IMPLICIT NONE
    REAL*8,ALLOCATABLE:: UTW(:,:),UTdW(:),LA(:,:),eigvec(:,:,:)
    REAL*8,ALLOCATABLE:: CM(:,:),RM(:,:),RMdum(:,:),Gdum1(:,:),Gdum2(:,:)
    INTEGER dumi,i,j,k,kk,SegShift
    REAL*8 dum,resparam

    CALL RANDOM_SEED()
    CALL Init()
    allocate(D(Nsvd),Dout(Nsvd,lambdanum),M(Msvd,lambdanum))
    allocate(V(Msvd,minMNsvd),W(Msvd),U(Nsvd,Msvd))
    allocate(UTdW(Msvd),LA(minMNsvd,lambdanum))
    if(NRseis>0)allocate(normdat(NSTAcomp))
    if(NRgps>0)allocate(sigmaGPS(3,NRgps))

    write(*,*)'Reading SVD.dat ...'
    open(111,form='unformatted',FILE='SVD.dat')
    read(111)D
    if(NRseis>0)read(111)normdat
    if(NRgps>0)then
      read(111)normdatGPS
      read(111)sigmaGPS
    endif
    read(111)W
    V=0.d0
    do i=1,minMNsvd
      read(111)V(1:Msvd,i)
    enddo
    read(111)dumi
    if(dumi.ne.0)then
      write(*,*)'Error! Rerun SlipInvSVD1!'
      stop
    else
      write(*,*)'  (check ok)'
    endif

    write(*,*)'Writing 10 leading eigenvectors ...'
    open(151,FILE='eigenvectors1d.dat')
    do k=1,10
      do kk=1,NSeg
        allocate(eigvec(Ssvd,NL(kk),NW(kk)))
        SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        eigvec(:,:,:)=RESHAPE(v(SegShift+1:SegShift+NW(kk)*NL(kk)*Ssvd,k),(/Ssvd, NL(kk), NW(kk)/))
        do i=1,Ssvd
          write(151,'(1000E13.5)')(sum(eigvec(i,j,1:NW(kk))),j=1,NL(kk))
        enddo
        write(151,*)
        write(151,*)
        deallocate(eigvec)
      enddo
    enddo
    close(151)
    open(151,FILE='eigenvectors.dat')
    do i=1,Msvd
      write(151,'(1000E13.5)')(v(i,k),k=1,10)
    enddo
    close(151)

    read(111)U
    read(111)dumi
    if(dumi.ne.0)then
      write(*,*)'Error! Rerun SlipInvSVD1!'
      stop
    else
      write(*,*)'  (check ok)'
    endif
    close(111)

    if(lambdanum>1)then
      write(*,*)'Calculating ',lambdanum," m_tilde's ..."
    else
      write(*,*)'Calculating m_tilde ...'
    endif

    maxw=maxval(W)

    if(abs(smoothkoef)>0.d0)then
      do i=1,Msvd
          U(:,i)=U(:,i)/W(i)
      enddo
    else
      do i=1,Msvd
        if(W(i)>maxw/1.d10)then
          U(:,i)=U(:,i)/W(i)
        else
          U(:,i)=0.d0
        endif
      enddo
    endif

    if(abs(smoothkoef)>0.d0)then
      write(*,*)'  (using all eigenvectors due to the applied smoothing)'
!	    if(minSVDchoice==1)then
!        write(*,*)'Error! MinSVD cannot be used when smoothing is applied...'
!        stop
!      endif
!      lambdafrom=Msvd
!      lambdato=Msvd
      lambdafrom=minMNsvd
      lambdato=minMNsvd
    elseif(eigsumchoice==1)then
      lambdafrom=0
      do i=1,minMNsvd
        if(W(i)>maxw)stop 'Error! Singular values are not in descending order!'
        if(W(i)>=maxw/lambdalim)then
          lambdafrom=lambdafrom+1
        endif
      enddo
      lambdato=lambdafrom
      write(*,*)'  (no. of eigenvectors: ', lambdafrom,')'
      write(*,*)'  (',dble(lambdafrom)/dble(Msvd)*100.d0,'% of eigenvectors used)'
    endif

#ifdef MKL
    call dgemv('T',Nsvd,Msvd,1.d0,U,Nsvd,D,1,0.d0,UTdW,1)
    deallocate(U)
    LA=0.d0
    do i=1,lambdanum
      LA(1:lambdafrom+i-1,i)=1.d0
    enddo
    do i=1,lambdanum
      LA(1:minMNsvd,i)=LA(1:minMNsvd,i)*UTdW(1:minMNsvd)
    enddo
    call dgemm('N','N',Msvd,lambdanum,minMNsvd,1.d0,V,Msvd,LA,minMNsvd,0.d0,M,Msvd)
#else
    UTdW=matmul(D,U)
    deallocate(U)
    LA=0.d0
    do i=1,lambdanum
      LA(1:lambdafrom+i-1,i)=1.d0
    enddo
    do i=1,lambdanum
      LA(1:minMNsvd,i)=LA(1:minMNsvd,i)*UTdW(1:minMNsvd)
    enddo
    M=matmul(V,LA)
#endif

!posterior covariance matrix CM   See (3.60) in Tarantola
    if(abs(smoothkoef)>0.d0)then
      write(*,*)'  (evaluating posterior covariance matrix CM...)'
      allocate(CM(Msvd,Msvd),RMdum(Msvd,Msvd),RM(Msvd,Msvd))
      do i=1,Msvd
        V(1:Msvd,i)=V(1:Msvd,i)/W(i)
      enddo
#ifdef MKL
      call dgemm('N','T',Msvd,Msvd,Msvd,1.d0,V,Msvd,V,Msvd,0.d0,CM,Msvd)
#else
      allocate(VT(minMNsvd,Msvd))
      VT=transpose(V)
      CM=matmul(V,VT)
      deallocate(VT)
#endif
      open(111,FILE='CM-stddev.dat')
      do i=1,Msvd
        write(111,*)sqrt(CM(i,i))
      enddo
      close(111)
      open(111,form='unformatted',FILE='CM.dat')
      write(111)CM
      close(111)
      open(111,FILE='mtildeslip2D-stddev.dat')
!      write(*,*)sqrt(sum(CM))*dt/dble(Msvd)   !standard deviation of the slip
      do kk=1,NSeg
        SegShift=sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        do j=1,NW(kk)
          dumi=SegShift+(j-1)*NL(kk)*Ssvd
          write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL(kk)),sqrt(sum(CM(dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd,dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd)))*dt  !støední odchylka skluzu
        enddo
        dumi=SegShift+(NW(kk)-1)*NL(kk)*Ssvd
        write(111,'(1000E13.5)')(sqrt(sum(CM(dumi+(i-1)*Ssvd+1:dumi+i*Ssvd,dumi+(i-1)*Ssvd+1:dumi+i*Ssvd)))*dt,i=1,NL(kk)),sqrt(sum(CM(dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd,dumi+(NL(kk)-1)*Ssvd+1:dumi+NL(kk)*Ssvd)))*dt  !støední odchylka skluzu
        close(111)
      enddo
    endif

    deallocate(V,LA)
    allocate(G(Nsvd,Msvd))
    open(111,form='unformatted',FILE='G.dat')
    do i=1,Msvd
      read(111)G(1:Nsvd,i)
    enddo
    close(111)
#ifdef MKL
    call dgemm('N','N',Nsvd,lambdanum,Msvd,1.d0,G,Nsvd,M,Msvd,0.d0,Dout,Nsvd)
#else
    Dout=matmul(G,M)
#endif

!resolution matrix RM     ! See (3.63) in Tarantola
    if(abs(smoothkoef)>0.d0)then
      write(*,*)'  (evaluating resolution matrix RM)'
      RM=0.d0;do i=1,Msvd;RM(i,i)=1.d0;enddo
#ifdef MKL
      call dgemm('T','N',Msvd,Msvd,Msvd,1.d0,G(Nseis+Ngps+2,1),Nsvd,G(Nseis+Ngps+2,1),Nsvd,0.d0,RMdum,Msvd)
      call dgemm('N','N',Msvd,Msvd,Msvd,-1.d0,CM,Msvd,RMdum,Msvd,1.d0,RM,Msvd)
#else
      allocate(Gdum1(Nsmooth,Msvd),Gdum2(Msvd,Nsmooth))
      Gdum1=G(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,1:Msvd)
      Gdum2=transpose(Gdum1)
      RMdum=matmul(Gdum2,Gdum1)
      deallocate(Gdum1,Gdum2)
      allocate(Gdum1(Msvd,Msvd))
      Gdum1=matmul(CM,RMdum)
      RM=RM-Gdum1
      deallocate(Gdum1)
#endif
      resparam=0.d0;do i=1,Msvd;resparam=resparam+RM(i,i);enddo
      write(*,*)'  (number of resolved parameters: ',resparam/dble(Msvd),' %)'
      open(111,FILE='RM.dat')
      do i=1,Msvd
        write(111,*)RM(i,i)
      enddo
      close(111)
      deallocate(RM,RMdum,CM)
    endif

    call OutputModel()

    open(201,FILE='UTdW.dat')
    do i=1,Msvd
      write(201,'(3E13.5)')maxw/W(i),UTdW(i),W(i)
    enddo
    close(201)

    END
