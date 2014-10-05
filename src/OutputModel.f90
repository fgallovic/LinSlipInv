!   AUTHOR:  Frantisek Gallovic
    
    SUBROUTINE OutputModel()
    USE SISVDmodule
    IMPLICIT NONE
    REAL*8,ALLOCATABLE:: varred(:),Mtilde(:,:,:),VRgps(:)
    REAL*8 M0tilde
    INTEGER i,j,k,kk
    integer,dimension(8) :: date
    REAL*4 stalocGPS(2)

    allocate(varred(lambdanum),VRgps(lambdanum),Mtilde(Ssvd,NL,NW))

    open(201,FILE='mtilde.dat')
    do j=1,Msvd
      write(201,'(1000E13.5)')(M(j,i),i=1,lambdanum)
    enddo
    close(201)

    if(NRseis>0)then
      do i=1,NSTAcomp
        if(abs(smoothkoef)>0.d0)then
          D((i-1)*nT+1:i*nT)=D((i-1)*nT+1:i*nT)*normdat(i)*sum(1.d0/normdat)/sum(staweight(:,:)*dble(stainfo(:,:)))*smoothkoefGF
          Dout((i-1)*nT+1:i*nT,:)=Dout((i-1)*nT+1:i*nT,:)*normdat(i)*sum(1.d0/normdat)/sum(staweight(:,:)*dble(stainfo(:,:)))*smoothkoefGF
        else
          D((i-1)*nT+1:i*nT)=D((i-1)*nT+1:i*nT)*normdat(i)
          Dout((i-1)*nT+1:i*nT,:)=Dout((i-1)*nT+1:i*nT,:)*normdat(i)
        endif
      enddo
      open(296,FILE='rvseisnez.dat')
      do i=1,nT
        write(296,'(1000E13.5)')dt*(iT1-1+i-1),(D((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(296)

      do i=1,lambdanum
        varred(i)=1.d0-sum((D(1:Nseis)-Dout(1:Nseis,i))**2)/sum(D(1:Nseis)**2)
      enddo

      open(484,FILE='varred.dat')
      open(297,FILE='svseisnez.dat')
    endif

    if(NRgps>0)then
      open(485,FILE='varredGPS.dat')
      open(298,FILE='mtildeslip2D-sGPS.out')
      open(292,file='stations-GPS-data.dat')
    endif

    open(296,FILE='mtildeslip2D.dat')
    open(201,FILE='mtilde1d.dat')
    open(202,FILE='mtilde.gnuplot.dat')
    open(295,FILE='mtildemomentrate.dat')
    open(299,FILE='srcmod.dat')

    do k=1,lambdanum
      Mtilde=RESHAPE(M(:,k),(/Ssvd, NL, NW/))
      do i=1,Ssvd
        write(201,'(1000E13.5)')(sum(Mtilde(i,j,:)),j=1,NL)
      enddo
      write(201,*);write(201,*)

      do i=1,Ssvd
        do j=1,NW
          write(202,'(1000E13.5)')Mtilde(i,1:NL,j)
        enddo
        write(202,*);write(202,*)
      enddo
      if(lambdanum==1)then
        if(NRseis>0)write(*,*)'  (variance reduction: ',varred(k),')'
        M0tilde=0.d0
        do i=1,Ssvd
          M0tilde=M0tilde+sum(mu(:,:)*Mtilde(i,:,:))*elem*dt
        enddo
        write(*,*)'  (scalar moment ',M0tilde,')'
        write(*,*)'  (scalar moment discrepancy ',Mfix/M0tilde*100.d0,'%)'
        if(smoothkoef>0.d0)then
          write(*,*)'  (RMS of covariance constraint ',sqrt(sum(Dout(Nseis+Ngps+2:Nseis+Ngps+1+Nsmooth,1)**2)/dble(Nsmooth))/abs(smoothkoef),')'
        endif
        call date_and_time(VALUES=date)
        write(299,'(A)')'# -----------------------------------------------------------------------------------------------------------'
        write(299,'(A)')'# SIV Inversion Exercise : xx'
        write(299,'(A,I2,A,I2,A,I4)')'# Date : ',date(3),'.',date(2),'.',date(1)
        write(299,'(A)')'# Modeler : F. Gallovic'
        write(299,'(A)')'# Inversion Method : Linear multi time-window with k^-2 smoothing'
        write(299,'(A)')'# Ground-motion code : Axitra'
        write(299,'(A,F8.3,A,E13.5)')'# SourcePar1 Mw-Mo [Nm] :',2./3.*log10(M0tilde)-6.0333,',',M0tilde
        write(299,'(A,F8.3,A,F8.3)')'# SourcePar2 L-W [km] :',leng/1000.,',',widt/1000.
        write(299,'(A,F8.3,A,F8.3)')'# SourcePar3 Strike-Dip [degrees] :',strike,',',dip
        write(299,'(A)')'# Hypocenter X-Y-Z [km] : 0, 0,', -hypodepth/1000.
        write(299,'(A,E13.5)')'# Depth2Top Z2top [km] :',-minval(sourZgps(:,:))/1000.+dW/2./1000.
        write(299,'(A,I,A,I)')'# NumPoints Nx-Nz :',NL,',',NW
        write(299,'(A,I,A,F8.3)')'# NumTimeWn Nt-Dt :',Ssvd,',',DT
        write(299,'(A)')'# ElemSTF : iso-tri '
        write(299,'(A)')'# ------------------------------------------------------------------------------------------------------------'
        write(299,'(A)')'# X            Y            Z            TotalSlip    Rake         RupTime      SlipTW       SlipTW       ...'
        write(299,'(A)')'# km           km           km           m            deg          s            m            m            ...'
        write(299,'(A)')'# ------------------------------------------------------------------------------------------------------------'
        do i=1,NL
          do j=1,NW
            write(299,'(1000E13.5)')sourEgps(i,j)/1000.,sourNgps(i,j)/1000.,-sourZgps(i,j)/1000.,sum(Mtilde(:,i,j))*dt,rakeGPS(i,j),T0,(Mtilde(kk,i,j)*dt,kk=1,Ssvd)
          enddo
        enddo
      endif

      do j=1,NW
        write(296,'(1000E13.5)')(sum(Mtilde(:,i,j))*dt,i=1,NL),sum(Mtilde(:,NL,j))*dt
      enddo
      write(296,'(1000E13.5)')(sum(Mtilde(:,i,NW))*dt,i=1,NL),sum(Mtilde(:,NL,NW))*dt
      write(296,*);write(296,*)

      do j=1,Ssvd
        write(295,'(1000E13.5)')dt*(j-1),sum(Mtilde(j,:,:)*mu(:,:))*dL*dW
      enddo
      write(295,*);write(295,*)

      if(NRseis>0)then
        write(484,'(I5,3E13.5)')lambdafrom+k-1,maxw/W(lambdafrom+k-1),smoothkoef,varred(k)
        do i=1,nT
          write(297,'(1000E13.5)')dt*(iT1-1+i-1),(Dout((j-1)*nT+i,k),j=1,NSTAcomp)
        enddo
        write(297,*);write(297,*)
      endif

      if(NRgps>0)then
        open(232,FILE='stations-GPS.dat',action='read')
        do i=1,NRgps
          read(232,*)stalocGPS(1:2)
          D(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3)=D(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3)*sigmaGPS(1:3,i)
          Dout(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3,k)=Dout(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3,k)*sigmaGPS(1:3,i)
          write(298,'(5E13.5)')stalocGPS(1:2),Dout(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3,k)
          write(292,'(5E13.5)')stalocGPS(1:2),D(Nseis+(i-1)*3+1:Nseis+(i-1)*3+3)
        enddo
        close(232)
        VRgps(k)=1.d0-sum((D(Nseis+1:Nseis+Ngps)-Dout(Nseis+1:Nseis+Ngps,k))**2)/sum(D(Nseis+1:Nseis+Ngps)**2)
        write(485,'(I5,3E13.5)')lambdafrom+k-1,maxw/W(lambdafrom+k-1),smoothkoef,VRgps(k)
      endif

    enddo

    deallocate(varred,Mtilde)

    END
